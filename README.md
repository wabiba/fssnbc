# FssNBC: Privacy-Preserving Naïve Bayes Classification via Function Secret Sharing

## 📖 Project Introduction

This project implements the FssNBC protocol based on [frankw2/libfss](https://github.com/frankw2/libfss).

FssNBC is an online efficient privacy-preserving Naïve Bayes classification protocol that utilizes Function Secret Sharing (FSS) technology to achieve low communication complexity and high computational efficiency while ensuring security.

Its core idea is:

- Use FSS to generate and evaluate secure distributed point functions (DPF/DCF) to support secure addition and comparison operations in probability calculations.

- Achieve the correctness and privacy of classification results without exposing the original training data and input samples.

### 🚀 Quick Start

```bash

# Execute in the root directory of the project

make fss-test

./fss-test
