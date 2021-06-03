---
layout: single 
title: "Preparing a Structure for Molecular Simulation"
categories:
  - Tutorial
tags:
  - Biophysics 
header: 
  image: https://salilab.org/modeller/gifs/modeller.jpg
---

<style type="text/css">
  body{
  font-size: 14pt;
}
</style>

## Prerequisites 

Prior to starting this tutorial, make sure you have the following installed: 

1. [Modeller](https://salilab.org/modeller/download_installation.html)
2. [Python](https://www.python.org/downloads/)

## Tutorial 

> **Background**

It is paramount that protein structures (.pdb) are prepared properly prior to conducting molecular simulations. If the quality of the sturcture is poor, any analysis from the simulation data will likely be useless. 

The most well-known database for the three-dimensional structural data of proteins is the [Protein Data Bank](https://www.rcsb.org/). This data is typically obtained by experimental techniques such as X-ray crystallography and cryo-electron microscopy. These techniques are not perfect and consequently there will be residues in the protein for which coordinate data is missing. 

This tutorial will walk through how we can use Modeller to **fill in gaps and introduce mutations for a protein structure**. 

> **Introduction**

Imagine that we are oncologists, and we are presented with a neuroblastoma patient who has a novel mutation, R1275Q, in the ALK gene. To investigate whether this mutation is the driver for this patient's neuroblastoma, we would like to conduct a MD simulation. Namely, we want to conduct an MD simulation of the ALK kinase in the inactive conformation, with the R1275Q mutation, and analyze changes in the hydrogen bond network of the protein. 

The first step is acquiring a structure of the inactive conformation of the ALK kinase (with no inhibitors bound). A quick google search shows us that the following is the entry of interest: [3LCS](https://www.rcsb.org/structure/3lcs). Create a directory for this tutorial, and save the PDB file in it. 

> **Filling in Missing Residues** 

As expected, upon inspection of the PDB file, we find that there are residues for which the coordinates are missing. We can follow Modeller's [tutorial](https://salilab.org/modeller/wiki/Missing%20residues) to fill in these gaps. 

```python
x = 'hello, python world!'
print(x.split(' '))
```

