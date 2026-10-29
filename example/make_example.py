#!/usr/bin/env python3
"""Write a tiny PLINK dataset (ex.bed/.bim/.fam, 4 people, 5 SNPs) and a PGS Catalog-style
model (model.txt) for trying out secretscores."""
snps = [(1, "rs1", 1000, "A", "G"), (1, "rs2", 2000, "C", "T"), (2, "rs3", 500, "G", "A"),
        (2, "rs4", 900, "T", "C"), (3, "rs5", 1500, "A", "C")]
weights = {"rs1": 0.10, "rs2": -0.05, "rs3": 0.20, "rs4": 0.15, "rs5": -0.12}
# copies of allele A1 per SNP and person (None = missing)
a1_copies = [[2, 1, 0, 2], [1, 1, 2, 0], [0, 2, 1, None], [2, 0, 1, 2], [1, 2, 0, 1]]
plink_code = {2: 0b00, 1: 0b10, 0: 0b11, None: 0b01}   # PLINK 2-bit genotype coding

with open("ex.bed", "wb") as f:
    f.write(bytes([0x6c, 0x1b, 0x01]))                  # magic number, SNP-major
    for row in a1_copies:                               # 4 people fit in one byte
        byte = 0
        for j, g in enumerate(row):
            byte |= plink_code[g] << (2 * j)
        f.write(bytes([byte]))
with open("ex.bim", "w") as f:
    for chrom, name, pos, a1, a2 in snps:
        f.write(f"{chrom}\t{name}\t0\t{pos}\t{a1}\t{a2}\n")
with open("ex.fam", "w") as f:
    for i in range(1, 5):
        f.write(f"fam{i} ind{i} 0 0 {1 if i % 2 else 2} -9\n")
with open("model.txt", "w") as f:
    f.write("#pgs_id=PGS000000\n#genome_build=GRCh37\n")
    f.write("rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n")
    for chrom, name, pos, a1, a2 in snps:
        f.write(f"{name}\t{chrom}\t{pos}\t{a1}\t{a2}\t{weights[name]}\t{chrom}\t{pos}\n")
