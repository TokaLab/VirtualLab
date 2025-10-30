# TokaLab - Virtual Laboratory

**TokaLab** is an open-access, open-source GitHub project designed to support multiple objectives:

* **Collaboration and knowledge sharing** – Provide a common platform where researchers working on nuclear fusion and tokamak physics can collaborate, share algorithms, and exchange knowledge.  
* **Education and training** – Create an educational environment for students and researchers by offering computational tools with different levels of physical fidelity.  
* **Data and algorithm development** – Build a flexible framework for data generation, algorithm validation, and machine learning model training.

---

## 🔬 VirtualLab - Overview

**VirtualLab** is the core component of TokaLab.  
It enables users to build a custom virtual tokamak or use existing configurations, generate plasma scenarios, and develop and compute synthetic diagnostics.  

It is designed to be **modular**, **object-oriented**, and to follow a **multi-fidelity physics** approach — ranging from simple educational codes to more advanced, research-grade tools.

Currently, two main modules are implemented:
* **SimPLa** (*Simulated Plasma*) – Solves the Grad–Shafranov equation on a fixed boundary.  
* **SynDiag** (*Synthetic Diagnostics*) – Includes codes for various diagnostic simulations.  

**VirtualLab** is actively maintained in both **MATLAB** and **Python**.

You can also use the MATLAB version directly online:

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=TokaLab/VirtualLab)

---

## 📚 Documentation

See our [Wiki](https://github.com/TokaLab/VirtualLab/wiki) for setup instructions, examples, and detailed explanations.

---

## 🗂 Repository Structure

```plaintext
VirtualLab/
│
├── VirtualLab_MATLAB/           # MATLAB
│   ├── docs
│   ├── examples
│   ├── SimPla_MATLAB
│   ├── SynDiag_MATLAB
│   ├── Validation            
│   └── VirtualLab_init.m
│
├── VirtualLab_Python/           # Python
│   └── ...
│ 
├── Citations.md
├── Contributing.md
├── Contributors.md
├── License
└── README.md             
```

---
## 🤝 Contributing

TokaLab is open-access and open-source — we warmly welcome new contributors!
Please check the [Contributing](./Contributing.md) guidelines to learn how to get started.

---
## 📄 License

TokaLab is released under the BSD 3-Clause License.
See the [License](./License) file for full details.

---
## 👥 TokaLab Team
The team and list of contributors are available [here](./Contributors.md).  

--- 
## 🧪 Research Outputs

TokaLab has already been used in research activities such as inverse problem algorithm validation and machine learning model training and testing. Explore our [Research Outputs](https://tokalab.github.io/Publications/) for more details.

---
## 📬 Contact

For questions, suggestions, or collaborations:

📧 Email: [tokalab.fusion@gmail.com](mailto:tokalab.fusion@gmail.com)  
🌐 Website: [tokalab.github.io](https://tokalab.github.io/)  
💼 Social: [LinkedIn](https://www.linkedin.com/company/tokalab-fusion/)  

