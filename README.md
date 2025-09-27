# 📉 Does the Shape of Economic Recovery Matter?  
### An Alternative Unit Root Test with New Smooth Transition Model  

This repository accompanies the article:  

**Özcan, M. (2022). Does the shape of economic recovery matter? An alternative unit root test with new smooth transition model. _The Journal of Economic Asymmetries, 26_, e00256. https://doi.org/10.1016/j.jeca.2022.e00256**

---

## 📘 About the Study  

The COVID-19 pandemic triggered sharp disruptions in global economic activity and gave rise to debates over the **shape of recovery** (V-, U-, or pit-shaped). This study introduces a **new smooth transition model** and an **alternative unit root testing procedure** to evaluate whether major macroeconomic indicators follow stationary or non-stationary processes in the presence of pandemic-induced structural changes.  

Key contributions of the study include:  
- Development of a **pit-shaped smooth transition function** to capture the nonlinear adjustment between pre- and post-pandemic trends.  
- Proposal of new **unit root test statistics** derived from this transition model.  
- **Monte Carlo simulations** to evaluate size and power properties of the new test.  
- An empirical application using **monthly industrial production index (IPI)**, **consumer price index (CPI)**, and **unemployment rates** of the G8 countries (2013–2021).  

**Findings:**  
The results demonstrate that pandemic-related shocks produced distinct pit-shaped adjustments in industrial production and unemployment across G8 economies, while consumer prices showed more heterogeneous responses. The new test proves to be more powerful than conventional linear tests in detecting such nonlinear dynamics:contentReference[oaicite:1]{index=1}.  

---

## 🛠 Repository Contents  

This repository contains **R codes** and **data files** to replicate the simulations and empirical analyses:  

- `Application.r` → Empirical application of the new test to IPI, CPI, and unemployment data.  
- `Criticals.r` → Simulation code for generating critical values.  
- `Powers.r` → Monte Carlo experiments for empirical power analysis.  
- `Sizes.r` → Monte Carlo experiments for finite sample size evaluation.  
- `plots.r` → Code for reproducing the figures in the article (smooth transitions and fitted trends).  
- `ipi.txt` → Industrial Production Index (IPI) data.  
- `cpi.txt` → Consumer Price Index (CPI) data.  
- `unp.txt` → Unemployment rate data.  
- `data.xlsx` → Additional dataset (compiled series used in the empirical application).  

All R scripts are structured for direct execution and can be adapted for further research.  

---

## 📂 Data Sources  

- **OECD Short-Term Economic Indicators** → Monthly IPI, CPI, and unemployment series for G8 countries (2013–2021).  
- **Author’s compilation** → Structured datasets (`.txt` and `.xlsx` files) used for empirical implementation.  

---

## ✨ Citation  

If you use this repository, please cite the article as:  

Özcan, M. (2022). *Does the shape of economic recovery matter? An alternative unit root test with new smooth transition model.* _The Journal of Economic Asymmetries, 26_, e00256. https://doi.org/10.1016/j.jeca.2022.e00256  

---

## 📝 License  

This repository is provided for **academic and research purposes only**. Please acknowledge the author when using or adapting the provided materials.  
