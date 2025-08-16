# NeuroSpikeLFP Toolkit

👩‍💻 Author - Wen  
- Email: annewen@zju.edu.cn  

## 📖 Introduction

This MATLAB-based toolkit integrates **spike merging, spike+LFP analysis, and visualization**, extending the functionality of **Chronux** for neural data research.  

---

## 📦 Modules

### 1. `neuroCombine`
- For different runs recorded at the same recording site, we use “neuroCombine” software combined with SUA's average waveform ISI、autocorrelation、3-dimentional PCA reasults and merge similar SUAs.
-example:
![neuroCombine](https://github.com/user-attachments/assets/f9650bfd-a578-4600-b5f9-2897ca03752e)

### 2. `neuroAnalysis`
- Spike, LFP and multi-scale batch analysis (bilt on Chronux)
- Function：
  - Spike features   
  - PSD for LFP
  - STA & FTA
  - Spike-LFP Coherence & Granger Causality
  -example:
<img width="1514" height="840" alt="1755316440375" src="https://github.com/user-attachments/assets/67a39c0b-df45-470e-b8d3-793cf5d2f517" />


### 3. `neuroViewer`
- Visualization of analysis results
-gui:
<img width="1349" height="671" alt="1755316716558" src="https://github.com/user-attachments/assets/543daea2-f48e-4caf-876e-95f1154058ca" />

-example results:
![neuroviewer](https://github.com/user-attachments/assets/6e58494d-525b-4ff8-b65d-6ef1d9362ebb)


## Requirements

- MATLAB R2023b or higher version
- [Chronux toolbox](http://chronux.org/)  

After installing dependencies, u can use this toolkit.
