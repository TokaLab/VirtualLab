# VirtualLab: Integrated Simulation and Diagnostics Playground

Welcome to **VirtualLab**, the example and demonstration environment of **TokaLab**. VirtualLab provides hands-on examples that integrate **SimPla** (Simulated Plasma Repository) and **SynDiag** (Synthetic Diagnostics) to form a complete virtual experimental setup. Users can simulate plasma equilibria and compute corresponding diagnostic signals in a seamless workflow.

---

## 🔬 Overview

VirtualLab is designed to:

* Demonstrate realistic workflows using SimPla and SynDiag
* Offer ready-to-run examples of equilibrium reconstruction + synthetic diagnostics
* Provide templates for new experiments and simulation chains

It serves as an educational and research-oriented sandbox to prototype virtual tokamak scenarios.

---

## 🗂 Repository Structure

```plaintext
VirtualLab/
│
├── examples/             # End-to-end use cases
│   ├── case1_simple/     # Basic equilibrium + diagnostic example
│   ├── case2_diverted/   # Advanced configuration with diagnostics
│   └── ...
│
├── shared/               # Shared utilities and data
├── plots/                # Output figures and diagnostics
└── README.md             # This file
```

---

## 🔄 Workflow

A typical VirtualLab workflow consists of:

1. **Define separatrix and compute plasma equilibrium** using **SimPla**
2. **Generate synthetic diagnostics** using **SynDiag**
3. **Visualize and interpret** results

---

## 🤖 Getting Started

1. Ensure `SimPla` and `SynDiag` are installed or accessible
2. Navigate to a case study under `examples/`
3. Run the provided script in Python or MATLAB

**Python Example:**

```bash
cd examples/case1_simple/
python run_case.py
```

**MATLAB Example:**

```matlab
cd examples/case1_simple/
run('run_case.m')
```

---

## 🔍 Included Use Cases

* **Circular Plasma with Interferometry**
* **D-Shaped Plasma with Thomson Scattering**
* **Diverted Configuration with Magnetic Diagnostics**

Each case includes:

* Inputs: separatrix, machine geometry, diagnostic setup
* Scripts: equilibrium reconstruction and diagnostics
* Outputs: plots, signals, performance metrics

---

## 🎓 Educational Value

VirtualLab is especially useful in:

* Classroom environments
* Thesis projects
* Rapid prototyping for new diagnostic methods

---

## 🤝 Contributing

Want to add your own case or use VirtualLab in your research? We welcome contributions:

1. Fork this repository
2. Add your new example under `examples/`
3. Submit a pull request

---

## 📄 License

Tokalab is licensed under the BSD 3-Clause License.  
Please see the [License](./License) file for full details.

---

## 📬 Contact

**TokaLab Team**
Email: \[[your-email@example.com](mailto:your-email@example.com)]
Website: \[your-lab-or-group-website]

---

Happy experimenting in the VirtualLab!
