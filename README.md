# Optical mapping error prediction
Software for prediction of errors in optical mapping by Bionano Saphyr instrument. The software currently supports human genome only.

# Requirements
* Python 3.5 or newer
* NumPy / SciPy see https://www.scipy.org/install.html

# Example
* Minimal running example: python3 main.py all.bnx -g genome.cmap 
* Store fragment lengths distribution: python3 main.py all.bnx -g genome.cmap --store_dist

# Test Data
* For testing purposes, please see the webpage of Bionano Inc. https://bionanogenomics.com/library/datasets/
