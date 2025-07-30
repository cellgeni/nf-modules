import torch
from tqdm import tqdm
from pyro.distributions import MultivariateNormal
import numpy as np
import pandas as pd

import itertools
import pysnooper


def e_step(data, w, theta, sigma, N, K, print_training_progress):
    # data = torch.clamp(data, min=0)  # Ensure no negative values before log10
    class_probs = torch.ones(N, K)
    for k in tqdm(range(K), disable=not print_training_progress):
        dist = MultivariateNormal(theta[k], sigma)
        class_probs[:, k] = w[k] * torch.exp(dist.log_prob(data))
    class_prob_norm = class_probs.div(torch.sum(class_probs, dim=1, keepdim=True))
    # class_prob_norm[torch.isnan(class_prob_norm)] = 0
    return class_prob_norm


def barcodes_01_from_channels(barcodes_1234, C, R):
    K = barcodes_1234.shape[0]
    barcodes_01 = np.ones((K, C, R))
    for b in range(K):
        barcodes_01[b, :, :] = 1 * np.transpose(
            barcodes_1234[b, :].reshape(R, 1) == np.arange(1, C + 1)
        )
    return barcodes_01


# auxiliary functions required for decoding
def torch_format(numpy_array):
    D = numpy_array.shape[1] * numpy_array.shape[2]
    return (
        torch.tensor(numpy_array)
        .float()
        .transpose(1, 2)
        .reshape(numpy_array.shape[0], D)
    )


def mat_sqrt(A, D):
    try:
        U, S, V = torch.svd(A + 1e-3 * A.mean() * torch.rand(D, D))
    except Exception:
        U, S, V = torch.svd(A + 1e-2 * A.mean() * torch.rand(D, D))
    S_sqrt = torch.sqrt(S)
    return torch.mm(torch.mm(U, torch.diag(S_sqrt)), V.t())

    # def map_states(data, N, D, C, R, K, codes, batch_size=None, temperature=0):
    #     inferred_model = infer_discrete(
    #         poutine.replay(
    #             model_constrained_tensor,
    #             trace=poutine.trace(auto_guide).get_trace(data, N, D, C, R, K, codes),
    #         ),
    #         temperature=temperature,
    #         first_available_dim=-2,
    #     )
    #     trace = poutine.trace(inferred_model).get_trace(data, N, D, C, R, K, codes)
    #     return trace.nodes["z"]["value"]


@pysnooper.snoop()
def decode_postprocess(
    model_params_and_losses_path,
    estimate_bkg=True,
    print_training_progress=True,
    inf_ind=np.array([]),
    estimate_additional_barcodes=None,
    modify_bkg_prior=True,  # should be False when there is a lot of background
    # signal (eg pixel-wise decoding, a lot of noisy
    # boundary tiles)
    add_remaining_barcodes_prior=0.05,  # after model is estimated, infeasible
    # barcodes are used in the e-step with
    # given prior
):
    params = torch.load(model_params_and_losses_path)
    w_star = params["w_star"]
    sigma_star = params["sigma_star"]
    sigma_ro_star = params["sigma_ro_star"]
    sigma_ch_star = params["sigma_ch_star"]
    theta_star = params["theta_star"]
    codes_tr_consts_v_star = params["codes_tr_consts_v_star"]
    codes_tr_v_star = params["codes_tr_v_star"]
    losses = params["losses"]
    K = params["K"]
    C = params["C"]
    D = params["D"]
    N = params["N"]
    R = params["R"]
    codes = params["codes"]
    data_norm = params["data_norm"]

    # include background / any additional barcode in codebook
    if estimate_bkg:
        bkg_ind = codes.shape[0]
        codes = torch.cat((codes, torch.zeros(1, D)))
    else:
        bkg_ind = np.empty((0,), dtype=np.int32)

    if np.any(estimate_additional_barcodes is not None):
        inf_ind = codes.shape[0] + np.arange(estimate_additional_barcodes.shape[0])
        codes = torch.cat((codes, torch_format(estimate_additional_barcodes)))
    else:
        inf_ind = np.empty((0,), dtype=np.int32)

    # computing class probabilities with appropriate prior probabilities
    if w_star is None or not isinstance(w_star, torch.Tensor) or len(w_star.shape) != 1:
        raise ValueError(
            "w_star must be a 1-dimensional torch.Tensor and cannot be None."
        )
    if modify_bkg_prior and w_star.shape[0] > K:
        # making sure that the K barcode classes have higher prior in case there are more than K
        # classes
        w_star_mod = torch.cat(
            (w_star[0:K], w_star[0:K].min().repeat(w_star.shape[0] - K))
        )
        w_star_mod = w_star_mod / w_star_mod.sum()
    else:
        w_star_mod = w_star

    if add_remaining_barcodes_prior > 0:
        # all possible barcodes
        barcodes_1234 = np.array(
            [p for p in itertools.product(np.arange(1, C + 1), repeat=R)]
        )
        # all possible barcodes in the same format as codes
        codes_inf = np.array(
            torch_format(barcodes_01_from_channels(barcodes_1234, C, R)).cpu()
        )
        # add the bkg code at the beginning
        codes_inf = np.concatenate((np.zeros((1, D)), codes_inf))
        codes_cpu = codes.cpu()
        for b in range(codes_cpu.shape[0]):  # remove already existing codes
            r = np.array(codes_cpu[b, :], dtype=np.int32)
            if np.where(np.all(codes_inf == r, axis=1))[0].shape[0] != 0:
                i = np.reshape(np.where(np.all(codes_inf == r, axis=1)), (1,))[0]
                codes_inf = np.delete(codes_inf, i, axis=0)
        if not estimate_bkg:
            bkg_ind = codes_cpu.shape[0]
            inf_ind = np.append(
                inf_ind, codes_cpu.shape[0] + 1 + np.arange(codes_inf.shape[0])
            )
        else:
            inf_ind = np.append(
                inf_ind, codes_cpu.shape[0] + np.arange(codes_inf.shape[0])
            )
        print(inf_ind.shape)
        codes_inf = torch.tensor(codes_inf).float()
        alpha = 1 - add_remaining_barcodes_prior
        w_star_all = torch.cat(
            (
                alpha * w_star_mod,
                torch.tensor((1 - alpha) / codes_inf.shape[0]).repeat(
                    codes_inf.shape[0]
                ),
            )
        )
        class_probs_star = e_step(
            data_norm,
            w_star_all,
            torch.matmul(
                torch.cat((codes, codes_inf)) * codes_tr_v_star
                + codes_tr_consts_v_star.repeat(w_star_all.shape[0], 1),
                mat_sqrt(sigma_star, D),
            ),
            sigma_star,
            N,
            w_star_all.shape[0],
            print_training_progress,
        )
    else:
        class_probs_star = e_step(
            data_norm,
            w_star_mod,
            theta_star,
            sigma_star,
            N,
            codes.shape[0],
            print_training_progress,
        )

    print(inf_ind, inf_ind.shape, class_probs_star.shape)
    print(class_probs_star)

    # collapsing added barcodes
    class_probs_star_s = torch.cat(
        (
            torch.cat(
                (
                    class_probs_star[:, 0:K],
                    class_probs_star[:, bkg_ind].reshape((N, 1)),
                ),
                dim=1,
            ),
            torch.sum(class_probs_star[:, inf_ind], dim=1).reshape((N, 1)),
        ),
        dim=1,
    )
    inf_ind_s = inf_ind[0] if len(inf_ind) > 0 else None
    # adding another class if there are NaNs
    nan_spot_ind = torch.unique(
        (torch.isnan(class_probs_star_s)).nonzero(as_tuple=False)[:, 0]
    )
    if nan_spot_ind.shape[0] > 0:
        nan_class_ind = class_probs_star_s.shape[1]
        class_probs_star_s = torch.cat(
            (class_probs_star_s, torch.zeros((class_probs_star_s.shape[0], 1))), dim=1
        )
        class_probs_star_s[nan_spot_ind, :] = 0
        class_probs_star_s[nan_spot_ind, nan_class_ind] = 1
    else:
        nan_class_ind = np.empty((0,), dtype=np.int32)

    class_probs = class_probs_star_s.cpu().numpy()

    class_ind = {
        "genes": np.arange(K),
        "bkg": bkg_ind,
        "inf": inf_ind_s,
        "nan": nan_class_ind,
    }
    torch_params = {
        "w_star": w_star_mod.cpu(),
        "sigma_star": sigma_star.cpu(),
        "sigma_ro_star": sigma_ro_star.cpu(),
        "sigma_ch_star": sigma_ch_star.cpu(),
        "theta_star": theta_star.cpu(),
        "codes_tr_consts_v_star": codes_tr_consts_v_star.cpu(),
        "codes_tr_v_star": codes_tr_v_star.cpu(),
        "losses": losses,
    }

    return {
        # 'class_probs': A numpy array of shape (N, number_of_classes) where each row contains
        # the posterior probabilities for each class (e.g., genes, background, infeasible, NaN)
        # for a given spot.
        "class_probs": class_probs,
        "class_ind": class_ind,
        "params": torch_params,
    }


# function handling output of decoding
def decoding_output_to_dataframe(out, df_class_names, df_class_codes):
    val = out["class_probs"].max(axis=1)
    ind = out["class_probs"].argmax(axis=1)
    K = len(out["class_ind"]["genes"])
    decoded = ind + 1
    decoded[np.isin(ind, out["class_ind"]["inf"])] = K + 1  # inf class
    decoded[np.isin(ind, out["class_ind"]["bkg"])] = K + 2  # bkg class
    decoded[np.isin(ind, out["class_ind"]["nan"])] = K + 3  # NaN class
    decoded_spots_df = pd.DataFrame(columns=["Name", "Code", "Probability"])
    decoded_spots_df["Name"] = df_class_names[decoded - 1]
    decoded_spots_df["Code"] = df_class_codes[decoded - 1]
    decoded_spots_df["Probability"] = val
    return decoded_spots_df


# function creating a heatmap for plotting spatial patterns
def heatmap_pattern(decoded_df, name, grid=150, thr=0.7, plot_probs=True):
    if not "Probability" in decoded_df.columns:
        if not "Score" in decoded_df.columns:
            plot_probs = False
            x_coord = np.floor(
                decoded_df.X[(decoded_df.Name == name)].to_numpy(dtype=np.double) / grid
            ).astype(np.int32)
            y_coord = np.floor(
                decoded_df.Y[(decoded_df.Name == name)].to_numpy(dtype=np.double) / grid
            ).astype(np.int32)
        else:
            x_coord = np.floor(
                decoded_df.X[
                    (decoded_df.Name == name) & (decoded_df.Score > thr)
                ].to_numpy(dtype=np.double)
                / grid
            ).astype(np.int32)
            y_coord = np.floor(
                decoded_df.Y[
                    (decoded_df.Name == name) & (decoded_df.Score > thr)
                ].to_numpy(dtype=np.double)
                / grid
            ).astype(np.int32)
    else:
        x_coord = np.floor(
            decoded_df.X[
                (decoded_df.Name == name) & (decoded_df.Probability > thr)
            ].to_numpy(dtype=np.double)
            / grid
        ).astype(np.int32)
        y_coord = np.floor(
            decoded_df.Y[
                (decoded_df.Name == name) & (decoded_df.Probability > thr)
            ].to_numpy(dtype=np.double)
            / grid
        ).astype(np.int32)
    H = np.zeros(
        (
            int(np.ceil(decoded_df.Y.to_numpy(dtype=np.double).max() / grid)),
            int(np.ceil(decoded_df.X.to_numpy(dtype=np.double).max() / grid)),
        )
    )
    if plot_probs:
        if "Probability" in decoded_df.columns:
            prob = decoded_df.Probability[decoded_df.Name == name].to_numpy(
                dtype=np.double
            )
        elif "Score" in decoded_df.columns:
            prob = decoded_df.Score[decoded_df.Name == name].to_numpy(dtype=np.double)
        prob[prob < thr] = 0
        for i in range(len(x_coord)):
            H[y_coord[i], x_coord[i]] = H[y_coord[i], x_coord[i]] + prob[i]
    else:
        coords = np.concatenate(
            (y_coord.reshape((len(x_coord), 1)), x_coord.reshape((len(x_coord), 1))),
            axis=1,
        )
        coords_u, coords_c = np.unique(coords, axis=0, return_counts=True)
        H[coords_u[:, 0], coords_u[:, 1]] = coords_c
    return H
