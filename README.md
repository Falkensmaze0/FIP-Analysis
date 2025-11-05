FIP Biomarker Analysis
======================

This repository contains `FIP_ANALYSIS_SCRIPT.py`, a reproducible pipeline for reconstructing feline infectious peritonitis (FIP) biomarker measurements from summary statistics and generating exploratory outputs such as correlation heatmaps, scatter matrices, and Welch group comparisons.

Quick Start
-----------

1. **Install Python 3.8+** (via python.org, Homebrew, or pyenv).
2. **Create a virtual environment** (recommended):
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```
3. **Install dependencies**:
   ```bash
   pip install pandas numpy seaborn matplotlib scipy
   ```
4. **Place input data**: ensure `biomarkers.csv` is in the same directory as the script. The CSV must include the columns `parameter`, `group`, `sample_size`, `mean`, `se`, `min`, and `max` (case-insensitive).
5. **Run the analysis**:
   ```bash
   python FIP_ANALYSIS_SCRIPT.py
   ```

What the Script Does
--------------------

- Reconstructs individual-level samples from group-level summary statistics using a constrained normal sampling scheme.
- Performs Z-score standardization to make biomarkers comparable before correlation analysis.
- Produces, for each group (healthy, dry, wet):
  - Annotated Pearson correlation heatmaps (`fip_analysis_results/correlation_heatmap_<group>.png`).
  - Pairwise scatterplot matrices with KDE diagonals (`fip_analysis_results/pairwise_scatter_<group>.png`).
  - Correlation coefficient matrices (`fip_analysis_results/correlation_matrix_<group>.csv`).
- Conducts Welch's t-tests across groups and stores the comparisons (`fip_analysis_results/group_comparisons_welch.csv`).
- Flags potential outliers based on |Z| > 3 and writes them to `fip_analysis_results/outliers_detected.csv`.

Working Directory and Outputs
-----------------------------

- Outputs are written to `fip_analysis_results/` beside the script. If the location is not writable, the script will fall back to `~/fip_analysis_results/` and notify you.
- Existing contents with the same file names are overwritten on each run.
- To start fresh, delete the output directory before rerunning.

Customising the Analysis
------------------------

- **Changing outlier sensitivity**: adjust the threshold in `detect_outliers` to a different Z-score cutoff.
- **Filtering biomarkers**: subset the DataFrame after loading `biomarkers.csv` to focus on specific markers.
- **Adapting for new groups**: extend the `valid_groups` list and ensure the CSV contains the new categories.
- **Saving plots elsewhere**: edit `OUTPUT_DIR` near the top of the script to a custom path.

Common Issues
-------------

- **Missing columns**: the script validates the header and raises a descriptive error if required fields are absent or misspelled.
- **Non-numeric values**: rows with non-numeric statistics are coerced to `NaN`. Inspect the CSV if reconstruction fails.
- **Small sample sizes**: Welch's t-test results may be unstable with n < 5; interpret p-values cautiously.

Reproducing Results
-------------------

The script seeds NumPy's random generator (`np.random.seed(42)`) so reconstructed samples and downstream outputs remain identical across reruns, provided the input CSV is unchanged.

Support
-------

For questions or to report issues, open a GitHub issue or contact the project maintainer.

تحلیل Biomarkerهای FIP
======================

این مخزن شامل `FIP_ANALYSIS_SCRIPT.py` است، یک pipeline قابل تکرار برای بازسازی اندازه‌گیری‌های FIP biomarker از آمار خلاصه و تولید خروجی‌های اکتشافی مانند correlation heatmaps، scatter matrices و Welch group comparisons.

شروع سریع
---------

1. **نصب Python 3.8+** (از طریق python.org، Homebrew یا pyenv).
2. **ایجاد محیط مجازی** (پیشنهاد می‌شود):
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```
3. **نصب وابستگی‌ها**:
   ```bash
   pip install pandas numpy seaborn matplotlib scipy
   ```
4. **قرار دادن داده‌های ورودی**: مطمئن شوید که `biomarkers.csv` در همان پوشه‌ای که اسکریپت قرار دارد موجود است. فایل CSV باید شامل ستون‌های `parameter`, `group`, `sample_size`, `mean`, `se`, `min` و `max` باشد (حساس به حروف نباشد).
5. **اجرای تحلیل**:
   ```bash
   python FIP_ANALYSIS_SCRIPT.py
   ```

کاری که اسکریپت انجام می‌دهد
----------------------------

- نمونه‌های سطح فردی را از آمار خلاصه سطح گروه با استفاده از یک طرح نمونه‌گیری نرمال محدود بازسازی می‌کند.
- استانداردسازی Z-score را برای قابل مقایسه کردن biomarkers پیش از تحلیل همبستگی انجام می‌دهد.
- برای هر گروه (healthy، dry، wet) تولید می‌کند:
  - correlation heatmaps پیرسون همراه با یادداشت (`fip_analysis_results/correlation_heatmap_<group>.png`).
  - scatterplot matrices جفتی با KDE diagonals (`fip_analysis_results/pairwise_scatter_<group>.png`).
  - correlation coefficient matrices (`fip_analysis_results/correlation_matrix_<group>.csv`).
- Welch's t-tests را بین گروه‌ها انجام می‌دهد و نتایج را در `fip_analysis_results/group_comparisons_welch.csv` ذخیره می‌کند.
- بر اساس |Z| > 3 outlierهای احتمالی را علامت‌گذاری کرده و آن‌ها را در `fip_analysis_results/outliers_detected.csv` می‌نویسد.

پوشه کاری و خروجی‌ها
--------------------

- خروجی‌ها در `fip_analysis_results/` کنار اسکریپت ذخیره می‌شوند. اگر این مکان قابل نوشتن نباشد، اسکریپت به `~/fip_analysis_results/` برمی‌گردد و به شما اطلاع می‌دهد.
- محتواهایی که نام فایل یکسان دارند در هر اجرا بازنویسی می‌شوند.
- برای شروع از ابتدا، پیش از اجرای دوباره، پوشه خروجی را حذف کنید.

سفارشی‌سازی تحلیل
-----------------

- **تغییر حساسیت outlier**: آستانه را در `detect_outliers` به یک cutoff Z-score متفاوت تنظیم کنید.
- **فیلتر کردن biomarkers**: پس از بارگذاری `biomarkers.csv` DataFrame را برای تمرکز بر نشانگرهای خاص زیرمجموعه کنید.
- **انطباق برای گروه‌های جدید**: لیست `valid_groups` را گسترش دهید و مطمئن شوید CSV شامل دسته‌های جدید است.
- **ذخیره نمودارها در جای دیگر**: `OUTPUT_DIR` را در ابتدای اسکریپت به مسیر دلخواه تغییر دهید.

مشکلات رایج
-----------

- **ستون‌های مفقود**: اسکریپت هدر را اعتبارسنجی می‌کند و اگر ستون‌های موردنیاز وجود نداشته باشند یا اشتباه تایپ شده باشند، خطای توصیفی می‌دهد.
- **مقادیر غیرعددی**: ردیف‌هایی با آمار غیرعددی به `NaN` تبدیل می‌شوند؛ اگر بازسازی شکست خورد فایل CSV را بررسی کنید.
- **حجم نمونه کوچک**: نتایج Welch's t-test ممکن است با n < 5 ناپایدار باشند؛ مقدار p را با احتیاط تفسیر کنید.

بازتولید نتایج
--------------

اسکریپت مولد تصادفی NumPy (`np.random.seed(42)`) را مقداردهی اولیه می‌کند، بنابراین نمونه‌های بازسازی‌شده و خروجی‌های پایین‌دستی تا زمانی که فایل ورودی CSV تغییری نکند یکسان باقی می‌مانند.

پشتیبانی
--------

برای پرسش‌ها یا گزارش مشکلات، یک GitHub issue باز کنید یا با نگه‌دارنده پروژه تماس بگیرید.
