# CCDS Codon Usage Explorer

Interactive GitHub Pages site for exploring codon usage patterns across all human CCDS genes.

```
wget https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_nucleotide.20221027.fna.gz
gunzip CCDS_nucleotide.20221027.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.20221027.txt
wget https://raw.githubusercontent.com/animesh/scripts/c6fa60a89eb4c3e835c6b997228eb907aae8fc9f/codonusage.pl
perl codonusage.pl CCDS_nucleotide.20221027.fna CCDS.20221027.txt 2>0
```


**Features**
- **PCA Explorer** – WebGL scatter plot of ~35 k genes coloured by cluster; click any gene to see its full codon-usage bar chart
- **Codon Heatmap** – cluster-centroid × codon frequency heatmap (all 64 codons, sorted by amino acid)
- **Gene Table** – searchable/filterable DataTables view; filter by cluster, chromosome, strand


## Repository layout

```
├── .github/
│   └── workflows/
│       └── deploy.yml       ← GitHub Actions: build & deploy
├── data/
│   └── CCDS_codon_usage.tsv ← (generated with perl codonusage.pl CCDS_nucleotide.20221027.fna CCDS.20221027.txt, codonusage.pl is available at https://raw.githubusercontent.com/animesh/scripts/c6fa60a89eb4c3e835c6b997228eb907aae8fc9f/codonusage.pl)
├── scripts/
│   └── process.py           ← PCA + KMeans → docs/data.json + docs/codons.json (generated with python scripts\process.py)
├── index.html               ← single-page app
└── README.md
```

The `docs/` directory is generated automatically by GitHub Actions – **do not commit it manually**.


### Add the TSV data file

```bash
cp /path/to/CCDS_nucleotide_20221027_fna_CCDS_20221027_txt_aa.txt \
   data/CCDS_codon_usage.tsv
```

> **Note:** The file is ~16 MB – well within GitHub's 100 MB limit... if it ever gets larger
> If you want to keep it out of the repo history, add `data/*.tsv` to `.gitignore`
> and store it in [Git LFS](https://git-lfs.github.com/) or GitHub Releases instead.

### Pages

https://animesh.github.io/codon-explorer/

---

## Customising the clustering

Edit the top of `scripts/process.py`:

```python
N_CLUSTERS   = 12    # number of KMeans clusters
RANDOM_STATE = 42    # reproducibility seed
```

Re-push to trigger a rebuild.

---

## Local preview (optional)

```bash
# Install deps
pip install pandas numpy scikit-learn scipy

# Generate data
python scripts/process.py

# Serve locally (Python built-in server)
cd docs
python -m http.server 8080
# Open http://localhost:8080
```

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| Actions workflow fails with `ModuleNotFoundError` | Check `pip install` step in `deploy.yml` |
| Page loads but scatter plot is empty | Open browser console; likely `data.json` fetch failed (path mismatch) |
| `data.json` is not found | Make sure `docs/` was generated and uploaded by Actions |
| Large `codons.json` loads slowly | Normal – it's ~9 MB; the status badge turns green when ready |
| `git push` rejects data file (too large) | Use [Git LFS](https://git-lfs.github.com/) for the `.tsv` |

---

## Citation

CCDS data https://ftp.ncbi.nlm.nih.gov/pub/CCDS/
