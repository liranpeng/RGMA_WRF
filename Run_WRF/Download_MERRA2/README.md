# MERRA2 Data Download Utility

This folder contains scripts and resources to download **MERRA2 reanalysis data** for use in modeling and analysis workflows.

---

## 📂 Contents

| File | Description |
|------|-------------|
| `Download_MERRA.py` | Python script used to download MERRA2 data files from the NASA GES DISC or other data servers. |
| `missing_data.txt` | A plain-text list specifying which MERRA2 files are missing and need to be downloaded. Each line typically contains a filename or date to process. |
| `MERRA.name` | A file summarizing the MERRA2 variable names, metadata, or naming conventions used in the dataset. |

---

## ⚙️ How It Works

1. **File List Input**

   The Download_MERRA.py script reads `missing_data.txt` to identify all files or dates that require downloading.

2. **Download Process**

   `Download_MERRA.py` connects to the data server, authenticates if needed, and downloads each specified file.

3. **Variables**

   The MERRA2 files downloaded include key atmospheric and aerosol variables. A summary of variable names and conventions can be found in `MERRA.name`.

---

## 🚀 Quick Start

To download missing data:

```bash
python Download_MERRA.py

