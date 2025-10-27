#!/usr/bin/env bash

set -euo pipefail

resource_type=$1
resource_path=$2

if [[ -z "$resource_type" || -z "$resource_path" ]]; then
    echo "Usage: $0 <resource_type> <resource_path>"
    exit 1
fi

imeta ls $resource_type "$resource_path" | \
    grep -v "AVUs" | \
    sed -e 's/^attribute: //' \
        -e 's/^value: //' \
        -e 's/^units: //' | \
    awk '
    NR%4==1 { attr=$0 }
    NR%4==2 { val=$0 }
    NR%4==3 { unit=$0 }
    NR%4==0 { printf "%s,%s,%s\n", attr, val, unit; attr=""; val=""; unit="" }
    END { if (attr != "") printf "%s,%s,%s\n", attr, val, unit }
    '