# Automated Multi-Gaussian Deconvolution for XRD Analysis

![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Version](https://img.shields.io/badge/Version-1.0-orange.svg)

## 📌 Overview
This repository contains the official MATLAB source code for the automated X-Ray Diffraction (XRD) analysis tool developed and presented in our research paper: 

> **"Development of an Automated Multi-Gaussian Deconvolution Algorithm for Quantitative Crystallinity Analysis from X-Ray Diffraction Patterns"**

This tool provides an end-to-end, automated computational pipeline to quantify the **Degree of Crystallinity ($X_c$)**, **Amorphous Fraction ($X_a$)**, and estimate **Crystallite Size ($D$)** using the Scherrer equation. It is particularly optimized for complex semi-crystalline materials (e.g., polymer matrices, starch, sago composites) where broad amorphous halos heavily overlap with sharp crystalline peaks.

---

## ✨ Key Features
*   **Asymmetric Least Squares (ALS) Baseline Correction:** Robust, parameter-free background removal without cutting off broad amorphous features.
*   **Automated Peak Detection:** Utilizes Savitzky-Golay filtering and local maxima detection to intelligently estimate initial parameters for optimization.
*   **Non-Linear Multi-Gaussian Deconvolution:** Employs the Levenberg-Marquardt algorithm for precise curve fitting, separating crystalline peaks from the amorphous halo.
*   **Analytical Area Integration:** Calculates phase areas ($A_c$, $A_a$) mathematically for high-precision crystallinity indexing (up to 5 decimal places).
*   **Automated Scherrer Equation:** Automatically filters crystalline phases and calculates crystallite size (in nm) based on the true analytical FWHM.
*   **Publication-Ready Output:** Generates high-quality, aesthetic subplots and automatically exports detailed numerical results to an Excel spreadsheet (`.xlsx`).

---

## 🚀 Getting Started

### Prerequisites
*   **MATLAB** (R2021a or newer recommended).
*   **Curve Fitting Toolbox** (Required for `lsqcurvefit` optimization).
*   **Signal Processing Toolbox** (Required for `findpeaks` and `sgolayfilt`).

### Installation
1.  Clone this repository to your local machine:
    ```bash
    git clone [https://github.com/YourUsername/XRD-MultiGaussian-Deconvolution.git](https://github.com/YourUsername/XRD-MultiGaussian-Deconvolution.git)
    ```
2.  Open MATLAB and navigate to the cloned directory.

### How to Use
1.  Prepare your raw XRD data in a two-column Excel file (`.xlsx`). 
    *   **Column 1:** $2\theta$ Angle (degrees).
    *   **Column 2:** Intensity (a.u.).
2.  Place your Excel file in the same directory as the MATLAB script.
3.  Open the main script (`Ultimate_MultiGaussian_XRD.m`).
4.  Change the `filename` variable in **Section 1** to match your Excel file name:
    ```matlab
    filename = 'YourData.xlsx'; 
    ```
5.  Click **Run**. The script will automatically process the data, display the results in the Command Window, plot the graphs, and save an Excel report (`Hasil_Lengkap_MultiGaussian_Scherrer.xlsx`).

*(Note: If no file is provided, the script will automatically generate a dummy dataset of a Sago+ZnO composite to demonstrate the algorithm's capabilities).*

---

## 📊 Output Visualization Examples

The algorithm generates a 2-panel visualization:
1.  **Top Panel:** Displays the raw data, ALS background, the cumulative Multi-Gaussian fit, and color-coded filled areas separating the amorphous ($X_a$) and crystalline ($X_c$) phases.
2.  **Bottom Panel:** Isolates the pure crystalline peaks and annotates each peak with its exact $2\theta$ position and estimated crystallite size ($D$ in nm) derived from the Scherrer equation.

*(You can add an image here later by uploading a screenshot of your MATLAB plot to your repository and linking it like this: `![Output Example](path_to_image.png)`)*

---

## 🔬 Scientific Background
The methodology implemented in this code addresses the subjectivity and reproducibility issues often found in manual XRD profile fitting. 

The algorithm classifies peaks based on their Full Width at Half Maximum (FWHM):
*   **Crystalline Phase:** $\text{FWHM} \le 2.5^\circ$
*   **Amorphous Phase:** $\text{FWHM} > 2.5^\circ$

Crystallinity ($X_c$) is calculated as:
$$X_c = \left( \frac{A_c}{A_c + A_a} \right) \times 100\%$$

Crystallite Size ($D$) is estimated using the Scherrer equation:
$$D = \frac{K \lambda}{\beta \cos \theta}$$
Where $K = 0.9$ (shape factor) and $\lambda = 0.15406$ nm (Cu-K$\alpha$ radiation).

---

## 📄 Citation
If you use this algorithm or code in your research, please cite our paper:

---

## 👨‍💻 Authors
*   **Junianto Sesa** (Universitas Papua)
*   **Bayu Harnadi Nasrul** (Universitas Hasanuddin)
*   **Andi Tessiwoja Tenri Ola** (Universitas Lambung Mangkurat)
*   **Pryandi M. Tabaika** (Universitas Sebelas September)
*   **Nurul Fajri Ramadhani Tang** (Universitas Papua)

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
