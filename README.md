# nf-modules

## Prerequisites
- nf-core tools:  `module load cellgen/nf-core/3.2.0` 

## Usage
### List all modules
This will show all all the modules available in the repository.
```bash
nf-core modules -g https://github.com/cellgeni/nf-modules.git list remote
```

### List all subworkflows
This will show all the subworkflows available in the repository.
```bash
nf-core subworkflows -g https://github.com/cellgeni/nf-modules.git list remote
```

### Install a subworkflow
This will install a subworkflow, along with the used `modules`, from the repository.
so, you shouldn't manually install `modules` from the repository unless is not used in a subworkflow.
```bash
nf-core subworkflows -g https://github.com/cellgeni/nf-modules.git install <subworkflow_name>
```

### Install a module
If the module is not used in a subworkflow, you can install it directly from the repository.
```bash
nf-core modules -g https://github.com/cellgeni/nf-modules.git install <module_name>
```

### List installed modules
This will show all the modules installed in your local repository.
```bash
nf-core modules -g https://github.com/cellgeni/nf-modules.git list local
```

### Remove a module
If you want to remove a module, you can do it with the following command:
```bash
nf-core modules -g https://github.com/cellgeni/nf-modules.git remove <module_name>d
```

### List all installed subworkflows
This will show all the subworkflows installed in your local repository.
```bash 
nf-core subworkflows -g https://github.com/cellgeni/nf-modules.git list local
```

### Remove a subworkflow
If you want to remove a subworkflow, you can do it with the following command:
```bash 
nf-core subworkflows -g https://github.com/cellgeni/nf-modules.git remove <subworkflow_name>
```