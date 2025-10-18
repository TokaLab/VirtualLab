# VirtualLab: Integrated Simulation and Diagnostics Playground

Welcome to **VirtualLab**, the example and demonstration environment of **TokaLab**. VirtualLab provides hands-on examples that integrate **SimPla** (Simulated Plasma Repository) and **SynDiag** (Synthetic Diagnostics) to form a complete virtual experimental setup. Users can simulate plasma equilibria and compute corresponding diagnostic signals in a seamless workflow.

---

## 🔬 Overview

VirtualLab is designed to:

* Demonstrate realistic workflows using SimPla and SynDiag
* Offer ready-to-run examples of equilibrium reconstruction + synthetic diagnostics
* Provide templates for new experiments and simulation chains

It serves as an educational and research-oriented sandbox to prototype virtual tokamak scenarios.
VirtualLab is especially useful in:

* Educational environments
* Thesis projects
* Digital twin simulations of tokamaks

---

## 🗂 Repository Structure

```plaintext
VirtualLab/
│
├── VirtualLab_MATLAB/           # MATLAB implementation (object-oriented)
│   ├── ...            
│
│── VirtualLab_MATLAB_edu/       # MATLAB function-oriented code for education
│   ├── ...
│
├── VirtualLab_Python/           # Python implementation (object-oriented)
│   └── ...
│ 
├── docs/                        # Coming soon!
│ 
├── License
└── README.md             
```

---

## 🔄 Workflow

A typical VirtualLab workflow consists of:

1. **Define separatrix and compute plasma equilibrium** using **SimPla**
2. **Generate synthetic diagnostics** using **SynDiag**
3. **Visualize and interpret** results

---

## 🤖 Getting Started

1. Ensure `SimPla` and `SynDiag` are downloaded in the proper folders VirtualLab_MATLAB/SimPla_MATLAB and VirtualLab_MATLAB/SynDiag_MATLAB (or VirtualLab_Python/SimPla_Python and VirtualLab_Python/SynDiag_Python)
2. Navigate to a case study under `examples/`
3. Run the provided script in Python or MATLAB

---

## 🤝 Contributing

We welcome contributions from the community. To contribute, please contact us.

---

## 📄 License

Tokalab is licensed under the BSD 3-Clause License.  
Please see the [License](./License) file for full details.

---

## 📬 Contact

For questions, suggestions, or collaborations:

**TokaLab Team**
Email: \[tokalab.fusion@gmail.com](mailto:tokalab.fusion@gmail.com)
Website: \[tokalab.github.io](https://tokalab.github.io/)


Happy experimenting in the VirtualLab!!





[![Open in MATLAB Online](https://matlab.mathworks.com/open/github/v1/badge.svg)](https://matlab.mathworks.com/open/github/v1?repo=TokaLab/VirtualLab)
