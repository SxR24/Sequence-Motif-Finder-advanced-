# 🔬 Advanced Sequence Motif Finder

## 📌 Overview

The **Advanced Sequence Motif Finder** is a powerful **bioinformatics tool** built using **Streamlit** that allows users to search for motifs in **DNA, RNA, or protein sequences**. It supports **exact and approximate matching (with mismatch tolerance)** and provides **interactive visualizations**, including:

- **Scrollable sequence view with highlighted motifs**
- **Motif match table**
- **Motif location visualization using Plotly**
- **Heatmap representation of motif intensity**
- **Circular genome visualization using DNA Features Viewer**
- **PDF summary report generation**

This tool is ideal for **researchers, bioinformaticians, and students** working with sequence analysis.

---

## 🚀 Installation Guide

To run this app locally, you need to install the required dependencies. Follow these steps:

### **1️⃣ Clone the Repository**

```sh
git clone https://github.com/SxR24/Sequence-Motif-Finder.git
cd Sequence-Motif-Finder
```

### **2️⃣ Set Up a Virtual Environment (Optional but Recommended)**

```sh
python -m venv venv
source venv/bin/activate  # On macOS/Linux
venv\Scripts\activate  # On Windows
```

### **3️⃣ Install Required Dependencies**

```sh
pip install -r requirements.txt
```

#### **Required Python Modules:**

The following modules are required:

- `streamlit`
- `re` (built-in)
- `pandas`
- `plotly`
- `matplotlib`
- `dna_features_viewer`
- `fuzzysearch`
- `fpdf`

You can install them manually if needed:

```sh
pip install streamlit pandas plotly matplotlib dna_features_viewer fuzzysearch fpdf
```

### **4️⃣ Run the App**

```sh
streamlit run app.py
```

The app will launch in your browser at `http://localhost:8501`.

---

## 📖 How It Works

### **🔹 User Inputs (Sidebar)**

- **Enter DNA/RNA/Protein sequence** in the text area.
- **Specify motifs** (comma-separated, supports regex patterns).
- **Set mismatch tolerance** (0-3 mismatches allowed).
- Click **🔍 Find Motifs** to analyze the sequence.

### **🔹 Outputs & Features**

#### **1️⃣ Scrollable Sequence View**

- The sequence is displayed in a **scrollable, highlighted format**, where matched motifs are color-coded.

#### **2️⃣ Motif Match Table**

- Displays all found motifs, their **start and end positions**, and the **matched sequence fragment**.

#### **3️⃣ Motif Location Plot (Plotly)**

- A **graphical representation** of motif positions along the sequence.

#### **4️⃣ Heatmap Visualization**

- Uses **color intensity** to show motif occurrence frequency.

#### **5️⃣ Circular Genome Visualization**

- Displays **graphical motifs** on a circular genome using `dna_features_viewer`.

#### **6️⃣ PDF Summary Report**

- Generates a **detailed report** of motif matches, sequence length, and analysis.
- Click **📥 Download PDF Summary** to save the report.

---

## 🎯 Example Usage

**Input Example:**

```
Sequence: AGCTAGCTAGCTGATCGTAGCTAGCTA
Motifs: AGCT, GAT
Mismatch Tolerance: 1
```

**Output:**

- Motif `AGCT` found at multiple positions.
- `GAT` detected with allowed mismatches.
- Sequence highlights, plots, and heatmap generated.

---

## 🛠 Future Enhancements

- 🧬 **Support for protein sequences with amino acid motifs**
- 📊 **Export results to CSV/Excel**
- ⚙️ **Adjustable visualization parameters**

---

## 🤝 Contributing

Contributions are welcome! Feel free to fork the repository and submit a **pull request**.

---

## 📜 License

This project is licensed under the **MIT License**.

---

## ⭐ Acknowledgments

- **Streamlit** for the web framework
- **Plotly** for interactive graphs
- **DNA Features Viewer** for genome visualization

---

## 📧 Contact

For questions or feedback, feel free to reach out!

📌 **GitHub:** [SxR24](https://github.com/SxR24)
📌 **Email:** [sohilananth109@gmail.com](mailto:sohilananth109@gmail.com)

