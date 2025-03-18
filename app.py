import streamlit as st
import re
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter
from fuzzysearch import find_near_matches
import matplotlib.pyplot as plt
import seaborn as sns
from dna_features_viewer import GraphicFeature, GraphicRecord
from io import BytesIO
from fpdf import FPDF
from Bio import SeqIO, Seq
import io
import numpy as np
import logomaker  # for sequence logo generation
from scipy.stats import binomtest

# ---------------------------
# Helper Functions

def canonical_kmer(kmer):
    """Return the lexicographically smallest of the k-mer and its reverse complement."""
    rc = str(Seq.Seq(kmer).reverse_complement())
    return min(kmer, rc)

def hamming_distance(s1, s2):
    """Compute Hamming distance between two equal-length strings."""
    if len(s1) != len(s2):
        return float('inf')
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# ---------------------------
# Streamlit Page Configuration
st.set_page_config(page_title="Sequence Motif Finder", layout="wide")
st.title('🔬 Advanced Sequence Motif Finder')

# Disclaimer
st.markdown("""
### **Disclaimer:**
This tool is under active development and may produce minor errors.  
Please verify outputs and cross-check with established tools when necessary.  
For suggestions or to report issues, please <a href="mailto:sohilananth109@gmail.com?subject=Feedback%20for%20Advanced%20Sequence%20Motif%20Finder">email me</a> or visit my <a href="https://github.com/SxR24">GitHub</a>.
""", unsafe_allow_html=True)

# ---------------------------
# Sidebar: Mode Selection & Input
mode = st.sidebar.radio("Select Search Mode", (
    "User Defined Motif", 
    "Basic De Novo Motif Discovery", 
    "Advanced De Novo Motif Discovery"
))

st.sidebar.header('⚙️ Input Sequence')
st.sidebar.subheader('Upload a FASTA file or paste the sequence below')
uploaded_file = st.sidebar.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

sequence = ""
if uploaded_file:
    fasta_text = uploaded_file.getvalue().decode("utf-8")
    fasta_io = io.StringIO(fasta_text)
    for record in SeqIO.parse(fasta_io, "fasta"):
        sequence += str(record.seq)
    st.sidebar.success(f"✅ Loaded sequence from {uploaded_file.name}")
else:
    sequence = st.sidebar.text_area('Enter DNA/RNA/Protein sequence:', height=200)

# Mode-specific sidebar inputs
if mode == "User Defined Motif":
    motif_input = st.sidebar.text_input('Enter motifs (comma-separated, regex supported):')
    mismatch_tolerance = st.sidebar.slider("Allowed mismatches:", 0, 3, 0)
elif mode == "Basic De Novo Motif Discovery":
    k = st.sidebar.number_input("Motif length (k)", min_value=4, max_value=12, value=6, step=1)
    top_n = st.sidebar.number_input("Number of top motifs to display", min_value=1, max_value=10, value=5, step=1)
    mismatch_tolerance = 0
    motif_input = ""
else:  # Advanced De Novo Motif Discovery
    min_k = st.sidebar.number_input("Minimum k-mer length", min_value=4, max_value=12, value=5, step=1)
    max_k = st.sidebar.number_input("Maximum k-mer length", min_value=min_k, max_value=12, value=7, step=1)
    p_threshold = st.sidebar.number_input("p-value threshold", min_value=0.0, max_value=1.0, value=0.05, step=0.01)
    cluster_threshold = st.sidebar.number_input("Clustering Hamming distance threshold", min_value=0, max_value=5, value=1, step=1)
    mismatch_tolerance = 0
    motif_input = ""

# Option: Reverse Complement Search
reverse_complement = st.sidebar.checkbox("Search in Reverse Complement")

# Remove whitespace
sequence = re.sub(r'\s+', '', sequence)

# Define colors for visualization
motif_colors = ["#FF6347", "#FFD700", "#32CD32", "#1E90FF", "#FF69B4"]

# Compute full sequence reverse complement if selected
if reverse_complement:
    rev_complement_seq = str(Seq.Seq(sequence).reverse_complement())
else:
    rev_complement_seq = None

# ---------------------------
# Main Processing: When "Find Motifs" button is clicked
if st.sidebar.button('🔍 Find Motifs'):
    if not sequence:
        st.error('⚠️ Please enter a sequence.')
    else:
        matches_dict = {}  # Will hold: motif -> {'matches': [(start, end, motif, strand), ...], 'color': color}
        if mode == "User Defined Motif":
            if not motif_input:
                st.error('⚠️ Please enter at least one motif.')
            else:
                motifs = [m.strip() for m in motif_input.split(',')]
                for i, motif in enumerate(motifs):
                    color = motif_colors[i % len(motif_colors)]
                    if mismatch_tolerance == 0:
                        orig_matches = [(m.start(), m.end(), motif, "+") for m in re.finditer(motif, sequence)]
                    else:
                        orig_matches = [(m.start, m.end, motif, "+") for m in find_near_matches(motif, sequence, max_l_dist=mismatch_tolerance)]
                    all_matches = orig_matches
                    if reverse_complement and rev_complement_seq is not None:
                        if mismatch_tolerance == 0:
                            rc_matches = []
                            for m in re.finditer(motif, rev_complement_seq):
                                orig_start = len(sequence) - m.end()
                                orig_end = len(sequence) - m.start()
                                rc_matches.append((orig_start, orig_end, motif, "-"))
                        else:
                            rc_matches = []
                            for m in find_near_matches(motif, rev_complement_seq, max_l_dist=mismatch_tolerance):
                                orig_start = len(sequence) - m.end
                                orig_end = len(sequence) - m.start
                                rc_matches.append((orig_start, orig_end, motif, "-"))
                        all_matches = all_matches + rc_matches
                    if all_matches:
                        matches_dict[motif] = {'matches': all_matches, 'color': color}
        elif mode == "Basic De Novo Motif Discovery":
            kmer_counts = {}
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if reverse_complement:
                    canon = canonical_kmer(kmer)
                else:
                    canon = kmer
                kmer_counts[canon] = kmer_counts.get(canon, 0) + 1
            sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
            motifs = [item[0] for item in sorted_kmers[:top_n]]
            st.markdown(f"**Basic De Novo Discovery:** Top {top_n} motifs of length {k}: {', '.join(motifs)}")
            for i, motif in enumerate(motifs):
                color = motif_colors[i % len(motif_colors)]
                pattern = motif
                orig_matches = [(m.start(), m.start()+k, motif, "+") for m in re.finditer(pattern, sequence)]
                all_matches = orig_matches
                if reverse_complement and rev_complement_seq is not None:
                    if motif != str(Seq.Seq(motif).reverse_complement()):
                        pattern_rc = f'(?=({motif}|{str(Seq.Seq(motif).reverse_complement())}))'
                        rc_matches = []
                        for m in re.finditer(pattern_rc, rev_complement_seq):
                            orig_start = len(sequence) - m.end()
                            orig_end = len(sequence) - m.start()
                            rc_matches.append((orig_start, orig_end, motif, "-"))
                        all_matches = all_matches + rc_matches
                if all_matches:
                    matches_dict[motif] = {'matches': all_matches, 'color': color}
        else:  # Advanced De Novo Motif Discovery
            # Build background model
            nuc_counts = Counter(sequence)
            total = len(sequence)
            bg = {nuc: nuc_counts[nuc]/total for nuc in nuc_counts}
            kmer_data = {}
            for k_val in range(int(min_k), int(max_k)+1):
                for i in range(len(sequence) - k_val + 1):
                    kmer = sequence[i:i+k_val]
                    canon = canonical_kmer(kmer) if reverse_complement else kmer
                    if canon not in kmer_data:
                        kmer_data[canon] = {"count": 0, "length": k_val}
                    kmer_data[canon]["count"] += 1
            significant_motifs = []
            for kmer, data in kmer_data.items():
                count = data["count"]
                length = data["length"]
                prob = np.prod([bg.get(nuc, 0) for nuc in kmer])
                trials = len(sequence) - length + 1
                expected = prob * trials
                result = binomtest(count, n=trials, p=prob, alternative='greater')
                p_val = result.pvalue
                if p_val < p_threshold:
                    significant_motifs.append((kmer, count, p_val))
            significant_motifs = sorted(significant_motifs, key=lambda x: x[2])
            st.markdown("**Advanced De Novo Discovery:**")
            st.markdown(f"Identified {len(significant_motifs)} candidate motifs from k={int(min_k)} to {int(max_k)} (p-value threshold={p_threshold}).")
            clusters = []
            used = set()
            for motif, count, p_val in significant_motifs:
                if motif in used:
                    continue
                cluster = [motif]
                used.add(motif)
                for other, count2, p_val2 in significant_motifs:
                    if other not in used and len(other) == len(motif):
                        if hamming_distance(motif, other) <= cluster_threshold:
                            cluster.append(other)
                            used.add(other)
                clusters.append(cluster)
            cluster_reps = []
            for cluster in clusters:
                rep = max(cluster, key=lambda m: kmer_data[m]["count"])
                cluster_reps.append(rep)
            motifs = cluster_reps
            st.markdown(f"**Final Representative Motifs:** {', '.join(motifs)}")
            for i, motif in enumerate(motifs):
                color = motif_colors[i % len(motif_colors)]
                pattern = motif
                orig_matches = [(m.start(), m.start()+len(motif), motif, "+") for m in re.finditer(pattern, sequence)]
                all_matches = orig_matches
                if reverse_complement and rev_complement_seq is not None:
                    if motif != str(Seq.Seq(motif).reverse_complement()):
                        pattern_rc = f'(?=({motif}|{str(Seq.Seq(motif).reverse_complement())}))'
                        rc_matches = []
                        for m in re.finditer(pattern_rc, rev_complement_seq):
                            orig_start = len(sequence) - m.end()
                            orig_end = len(sequence) - m.start()
                            rc_matches.append((orig_start, orig_end, motif, "-"))
                        all_matches = all_matches + rc_matches
                if all_matches:
                    matches_dict[motif] = {'matches': all_matches, 'color': color}

        # ---------------------------
        # Highlighted Sequence View
        positions = []
        for motif, data in matches_dict.items():
            color = data['color']
            for start, end, _, strand in data['matches']:
                positions.append((start, f"<span style='background-color:{color}; padding:2px; border-radius:4px;'>"))
                positions.append((end, "</span>"))
        positions.sort(key=lambda x: x[0])
        highlighted_sequence = ""
        prev_index = 0
        for index, tag in positions:
            highlighted_sequence += sequence[prev_index:index] + tag
            prev_index = index
        highlighted_sequence += sequence[prev_index:]
        st.subheader("📜 Scrollable Sequence View")
        st.markdown(f"""
            <div style='overflow-x: auto; white-space: nowrap; font-family: monospace; font-size: 16px;'>
                {highlighted_sequence}
            </div>
        """, unsafe_allow_html=True)
        st.markdown("""
**Scientific Explanation:**  
The highlighted sequence view visually maps the discovered motifs onto your sequence,  
allowing you to inspect their spatial distribution and potential functional relationships.
""")
        # ---------------------------
        # Results Table
        results = []
        for motif, data in matches_dict.items():
            for start, end, _, strand in data['matches']:
                if strand == "+":
                    matched_seq = sequence[start:end]
                else:
                    matched_seq = str(Seq.Seq(sequence[start:end]).reverse_complement())
                results.append([motif, start, end, matched_seq, strand])
        if results:
            df_results = pd.DataFrame(results, columns=['Motif', 'Start', 'End', 'Matched Sequence', 'Strand'])
            st.subheader('🔍 Motif Matches')
            st.dataframe(df_results)
            st.markdown("""
**Scientific Explanation:**  
The results table provides quantitative details (positions, strand, sequence) of each motif occurrence,  
serving as a basis for further statistical or biological interpretation.
""")
        # ---------------------------
        # PDF Report Generation
        def generate_summary_pdf():
            pdf = FPDF()
            pdf.set_auto_page_break(auto=True, margin=15)
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            pdf.cell(200, 10, "Sequence Motif Finder - Summary Report", ln=True, align='C')
            pdf.ln(10)
            pdf.cell(200, 10, f"Total Sequence Length: {len(sequence)}", ln=True)
            if mode == "User Defined Motif":
                pdf.cell(200, 10, f"Motifs Searched: {', '.join(motifs)}", ln=True)
            elif mode == "Basic De Novo Motif Discovery":
                pdf.cell(200, 10, f"De Novo Discovered Motifs: {', '.join(motifs)}", ln=True)
            else:
                pdf.cell(200, 10, f"Advanced De Novo Representative Motifs: {', '.join(motifs)}", ln=True)
            pdf.cell(200, 10, f"Total Matches Found: {sum(len(data['matches']) for data in matches_dict.values())}", ln=True)
            pdf.ln(10)
            pdf.cell(200, 10, "Motif Matches:", ln=True)
            pdf.ln(5)
            for motif, data in matches_dict.items():
                for start, end, _, strand in data['matches']:
                    if strand == "+":
                        seq_match = sequence[start:end]
                    else:
                        seq_match = str(Seq.Seq(sequence[start:end]).reverse_complement())
                    pdf.cell(200, 10, f"Motif: {motif}, Start: {start}, End: {end}, Strand: {strand}, Sequence: {seq_match}", ln=True)
            return pdf.output(dest='S').encode('latin1')
        
        st.subheader("📄 Download Summary Report")
        pdf_data = generate_summary_pdf()
        st.download_button("📥 Click here to Download PDF Summary", data=BytesIO(pdf_data), file_name="Motif_Summary_Report.pdf", mime="application/pdf")
        st.markdown("""
**Scientific Explanation:**  
The PDF report compiles all motif search parameters, results, and statistics into a permanent record  
that can be used for documentation, sharing, and further analysis.
""")
        # ---------------------------
        # 2D Visualization & Heatmap
        st.subheader('📈 Motif Visualization & Heatmap')
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=list(range(len(sequence))),
            y=[0]*len(sequence),
            mode='lines',
            line=dict(color='lightgrey', width=2),
            name='Sequence'
        ))
        heatmap_data = [0] * len(sequence)
        for motif, data in matches_dict.items():
            color = data['color']
            for start, end, _, _ in data['matches']:
                for i in range(start, end):
                    heatmap_data[i] += 1
                fig.add_trace(go.Scatter(
                    x=list(range(start, end)),
                    y=[0.5]*len(range(start, end)),
                    mode='markers+lines',
                    marker=dict(size=10, color=color, opacity=0.8),
                    line=dict(color=color, width=4),
                    name=f'{motif} (Pos {start}-{end})'
                ))
        fig.update_layout(
            title='🔬 Motif Locations in Sequence',
            xaxis_title='Sequence Position',
            yaxis_title='Motif Intensity',
            hovermode='x unified'
        )
        st.plotly_chart(fig, use_container_width=True)
        st.markdown("""
**Scientific Explanation:**  
This 2D plot visually represents where motifs occur along the sequence and the intensity (overlap) of occurrences,  
helping to pinpoint regions of potential regulatory importance.
""")
        # ---------------------------
        st.subheader("🔥 Motif Frequency Bar Chart")
        motif_counts = Counter([motif for motif, _ in matches_dict.items()])
        fig_bar = px.bar(x=list(motif_counts.keys()), y=list(motif_counts.values()),
                         labels={'x':'Motif', 'y':'Count'}, title='Motif Frequency')
        st.plotly_chart(fig_bar, use_container_width=True)
        st.markdown("""
**Scientific Explanation:**  
The bar chart provides a quantitative comparison of motif prevalence, offering insights into which motifs  
may play dominant regulatory or structural roles.
""")
        # ---------------------------
        st.subheader("📌 Motif Distribution Across Sequence")
        fig_line = px.line(x=list(range(len(sequence))), y=heatmap_data,
                           labels={'x':'Position', 'y':'Motif Occurrences'})
        st.plotly_chart(fig_line, use_container_width=True)
        st.markdown("""
**Scientific Explanation:**  
This line plot maps the overall distribution of motifs along the sequence, revealing hotspots  
of motif enrichment that might be associated with functional genomic regions.
""")
        # ---------------------------
        st.subheader("🧬 GC Content Visualization")
        gc_content = [(sequence[:i+1].count('G') + sequence[:i+1].count('C'))/(i+1) for i in range(len(sequence))]
        fig_gc = px.line(x=list(range(len(sequence))), y=gc_content,
                         labels={'x':'Position', 'y':'GC Content'})
        st.plotly_chart(fig_gc, use_container_width=True)
        st.markdown("""
**Scientific Explanation:**  
The GC content plot offers insights into the nucleotide composition of the sequence,  
which can affect its structural stability and protein-binding affinity.
""")
        # ---------------------------
        st.subheader("🔥 Heatmap Visualization")
        heatmap_fig = px.imshow([heatmap_data], color_continuous_scale='viridis', aspect='auto')
        st.plotly_chart(heatmap_fig, use_container_width=True)
        st.markdown("""
**Scientific Explanation:**  
The heatmap offers a visual overview of motif intensity along the sequence, helping to identify  
regions with significant motif enrichment or depletion.
""")
        # ---------------------------
        # Only show Motif Co-occurrence Heatmap and Sequence Logo Plot for User Defined Motif mode
        if mode == "User Defined Motif":
            st.subheader("🔥 Motif Co-occurrence Heatmap")
            co_matrix = pd.DataFrame([[1 if motif in sequence[start:end] else 0 for motif in motifs] 
                                      for _, start, end, _, _ in results], columns=motifs)
            fig_heat, ax = plt.subplots()
            sns.heatmap(co_matrix.corr(), annot=True, cmap='coolwarm', ax=ax)
            st.pyplot(fig_heat)
            st.markdown("""
**Scientific Explanation:**  
The co-occurrence heatmap displays correlations between motif occurrences,  
potentially indicating cooperative interactions or shared regulatory mechanisms.
""")
            st.subheader("🏆 Sequence Logo Plot")
            logo_figures = []
            for motif, data in matches_dict.items():
                seqs = []
                for start, end, _, strand in data['matches']:
                    if strand == "+":
                        seqs.append(sequence[start:end])
                    else:
                        seqs.append(str(Seq.Seq(sequence[start:end]).reverse_complement()))
                st.markdown(f"**Sequence Logo for motif: {motif}**")
                if len(seqs) > 1 and len(set(len(s) for s in seqs)) == 1:
                    count_mat = logomaker.alignment_to_matrix(sequences=seqs, to_type='counts')
                    fig_logo, ax = plt.subplots(figsize=(3, 2))
                    logomaker.Logo(count_mat, ax=ax)
                    ax.set_title(f"{motif}", fontsize=10)
                    ax.set_xticks([])  # Remove x-ticks
                    ax.set_yticks([])  # Remove y-ticks
                    logo_figures.append(fig_logo)
                else:
                    st.write(f"Not enough uniform matches for motif {motif} to generate a sequence logo.")
            if logo_figures:
                for i in range(0, len(logo_figures), 2):
                    cols = st.columns(2)
                    cols[0].pyplot(logo_figures[i])
                    if i + 1 < len(logo_figures):
                        cols[1].pyplot(logo_figures[i+1])
            st.markdown("""
**Scientific Explanation:**  
The sequence logo plot graphically represents the conservation and variability across aligned motif matches.  
Higher letter stacks at a given position indicate high information content, suggesting a functional importance.
""")
        # ---------------------------
        st.subheader("🌀 Interactive Circular Genome Visualization")
        fig_circular = go.Figure()
        theta_outer = np.linspace(0, 360, 361)
        r_outer = np.full_like(theta_outer, 1)
        fig_circular.add_trace(go.Scatterpolar(
            r=r_outer,
            theta=theta_outer,
            mode='lines',
            line=dict(color='black', width=2),
            name='Genome'
        ))
        for motif, data in matches_dict.items():
            color = data['color']
            for start, end, _, strand in data['matches']:
                start_angle = 360 * start / len(sequence)
                end_angle = 360 * end / len(sequence)
                theta_arc = np.linspace(start_angle, end_angle, 100)
                r_val = 1.2 if strand == "+" else 0.8
                r_arc = np.full_like(theta_arc, r_val)
                fig_circular.add_trace(go.Scatterpolar(
                    r=r_arc,
                    theta=theta_arc,
                    mode='lines',
                    line=dict(color=color, width=4),
                    name=f'{motif} ({strand})',
                    showlegend=True
                ))
        fig_circular.update_layout(
            polar=dict(
                radialaxis=dict(visible=True, range=[0, 1.5]),
                angularaxis=dict(direction="clockwise")
            ),
            showlegend=True,
            title="Interactive Circular Genome Visualization"
        )
        st.plotly_chart(fig_circular, use_container_width=True)
        st.markdown("""
**Scientific Explanation:**  
This polar plot maps motif occurrences onto a circular representation of the sequence,  
which is particularly valuable for studying the spatial organization of circular genomes (e.g., plasmids, mitochondrial DNA).
""")
        # ---------------------------
        st.subheader("📊 3D Motif Distribution Plot")
        motif_index = {motif: idx for idx, motif in enumerate(motifs)}
        x_vals, y_vals, z_vals, text_vals = [], [], [], []
        for motif, data in matches_dict.items():
            for start, end, _, strand in data['matches']:
                x_vals.append(start)
                y_vals.append(motif_index[motif])
                z_vals.append(end - start)
                text_vals.append(f"{motif} ({strand})<br>Pos: {start}-{end}")
        if x_vals:
            fig_3d = go.Figure(data=[go.Scatter3d(
                x=x_vals,
                y=y_vals,
                z=z_vals,
                mode='markers',
                marker=dict(
                    size=5,
                    color=z_vals,
                    colorscale='Viridis',
                    opacity=0.8
                ),
                text=text_vals,
                hoverinfo='text'
            )])
            fig_3d.update_layout(
                scene=dict(
                    xaxis_title='Start Position',
                    yaxis=dict(title='Motif', tickvals=list(motif_index.values()), ticktext=list(motif_index.keys())),
                    zaxis_title='Match Length'
                ),
                title="3D Motif Distribution Plot"
            )
            st.plotly_chart(fig_3d, use_container_width=True)
            st.markdown("""
**Scientific Explanation:**  
The 3D scatter plot offers a multidimensional view of motif occurrences:  
- **X-axis:** Start position of each motif occurrence.  
- **Y-axis:** Differentiates motif types.  
- **Z-axis:** Represents the match length, indicating potential structural variation.  
This comprehensive view aids in detecting complex spatial patterns in the data.
""")
