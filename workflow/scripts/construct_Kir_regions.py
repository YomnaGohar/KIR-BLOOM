#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 22:22:15 2025

@author: yomna
"""

# Input files
fai_file = "/home/yomna/hpc_project/kir/reference/refilt_ref_2/ref.fa.fai"
txt_file = "/home/yomna/sevenOfNine/cram_files/ERR3239430/chromsome_naming.txt"

# Read .fai (format: name length offset ...)
fai_lengths = {}
with open(fai_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            name, length = parts[0], int(parts[1])
            fai_lengths[length] = name
sequences_to_remove = {
    "GL000209.2", "NT_113949.2",
    "GL949746.1", "NW_003571054.1",
    "GL949747.2", "NW_003571055.2",
    "GL949748.2", "NW_003571056.2",
    "GL949749.2", "NW_003571057.2",
    "GL949750.2", "NW_003571058.2",
    "GL949751.2", "NW_003571059.2",
    "GL949752.1", "NW_003571060.1",
    "GL949753.2", "NW_003571061.2",
    "KI270882.1", "NT_187636.1",
    "KI270883.1", "NT_187637.1",
    "KI270884.1", "NT_187638.1",
    "KI270885.1", "NT_187639.1",
    "KI270886.1", "NT_187640.1",
    "KI270887.1", "NT_187641.1",
    "KI270888.1", "NT_187642.1",
    "KI270889.1", "NT_187643.1",
    "KI270930.1", "NT_187684.1",
    "KI270931.1", "NT_187685.1",
    "KI270932.1", "NT_187686.1",
    "KI270933.1", "NT_187687.1",
    "KI270938.1", "NT_187693.1",
    "KV575246.1", "NW_016107300.1",
    "KV575247.1", "NW_016107301.1",
    "KV575248.1", "NW_016107302.1",
    "KV575249.1", "NW_016107303.1",
    "KV575250.1", "NW_016107304.1",
    "KV575251.1", "NW_016107305.1",
    "KV575252.1", "NW_016107306.1",
    "KV575253.1", "NW_016107307.1",
    "KV575254.1", "NW_016107308.1",
    "KV575255.1", "NW_016107309.1",
    "KV575256.1", "NW_016107310.1",
    "KV575257.1", "NW_016107311.1",
    "KV575258.1", "NW_016107312.1",
    "KV575259.1", "NW_016107313.1",
    "KV575260.1", "NW_016107314.1",
    "KI270890.1", "NT_187644.1",
    "KI270891.1", "NT_187645.1",
    "KI270914.1", "NT_187668.1",
    "KI270915.1", "NT_187669.1",
    "KI270916.1", "NT_187670.1",
    "KI270917.1", "NT_187671.1",
    "KI270918.1", "NT_187672.1",
    "KI270919.1", "NT_187673.1",
    "KI270920.1", "NT_187674.1",
    "KI270921.1", "NT_187675.1",
    "KI270922.1", "NT_187676.1",
    "KI270923.1", "NT_187677.1",
    "KI270929.1", "NT_187683.1"
}
# Read txt (format: name length)
txt_lengths = {}
with open(txt_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            name, length = parts[0], int(parts[1])
            txt_lengths[length] = name
output_file="/home/yomna/sevenOfNine/cram_files/cram_files/resources/kir_regions.bed"
i = 0
with open(output_file, "w") as out:
    out.write(f"chr19\t52025634\t57084318\n")
    for length, fai_name in fai_lengths.items():
        if length in txt_lengths and fai_name in sequences_to_remove:
            txt_name = txt_lengths[length]
            out.write(f"{txt_name}\t0\t{length}\n")
            i += 1

print(f"Saved {i} alternate chr19 sequences to: {output_file}")
