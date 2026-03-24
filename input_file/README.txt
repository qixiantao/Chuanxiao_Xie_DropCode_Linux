# DEMO Data – Test Your Installation

This folder contains a **small example dataset** that you can use to verify that GEAtool has been installed correctly and all required tools are functioning.

## What's Inside

- **`sample1_1.fq.gz`** & **`sample1_2.fq.gz`**: Simulated paired‑end reads (gzipped).
- **`reference.fasta`**: A short reference sequence containing the target region.
- **`target.fasta`**: The target editing site used for extraction.
- **`barcode.xlsx`**: Example barcode mapping (two samples).
- **`library.xlsx`** (optional): Sample name to password mapping (if you want to test name replacement).
1. **Run the pipeline on the DEMO dataset**  
   The DEMO folder is already located under `input_file/DEMO/`. To process it:
   ```bash
   bash batch_run.sh --input ./input_file --output ./output_test
   ```
   (You can also use `--t`, `--ram`, `--q` to adjust resources.)

2. **Check the output**  
   After a successful run, you will find results in `output_test/DEMO/`:
   - `BAM/` folder with aligned BAM files and indices.
   - `GEAnalysis_result.xlsx` containing allele frequencies.

   If the Excel file is non‑empty and contains the expected top alleles, the pipeline is working correctly.

## What to Expect

The example data should produce:
- **~5000 reads** per sample after filtering.
- **Several allele variants** (e.g., A, C, G, T) with frequencies summing to 100%.
- **No errors** in the logs (you may see warnings about missing `library.xlsx` – that is normal).

## Important: Remove DEMO Before Your Real Analysis

**Before you start processing your own samples, please delete or move the DEMO folder to avoid mixing it with your real data:**

```bash
rm -rf ./input_file/DEMO
```

If you prefer to keep the demo for later reference, you can move it elsewhere:
```bash
mv ./input_file/DEMO ./DEMO_backup
```

## Troubleshooting

- If the pipeline fails, check the logs in `logs/` for details.
- If the Excel result is empty, ensure that `target.fasta` matches the reference and that the barcodes in `barcode.xlsx` are correct.
- For any other issues, refer to the main [README](../README.md) or open an issue on GitHub.

---

Happy testing! 🚀