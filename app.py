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

# Streamlit Page Configuration
st.set_page_config(page_title="Sequence Motif Finder", layout="wide")
st.title('🔬 Advanced Sequence Motif Finder')

# Disclaimer below the title
st.subheader("""
**Disclaimer:**
This tool is currently under development and may produce minor errors.  
Please verify the output carefully and cross-check with other tools if necessary.  
If you notice any logical errors or have suggestions for improvements, feel free to  
email me at sohilananth109@gmail.com or reach me at my github https://github.com/SxR24 .
""")

# Sidebar Inputs
st.sidebar.header('⚙️ Input Sequence')
st.sidebar.subheader('Upload FASTA file or paste the sequence in the box')

# File Uploader for FASTA Files
uploaded_file = st.sidebar.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

# Sequence Input Area
sequence = ""
if uploaded_file:
    fasta_text = uploaded_file.getvalue().decode("utf-8")
    fasta_io = io.StringIO(fasta_text)
    for record in SeqIO.parse(fasta_io, "fasta"):
        sequence += str(record.seq)
    st.sidebar.success(f"✅ Loaded sequence from {uploaded_file.name}")
else:
    sequence = st.sidebar.text_area('Enter DNA/RNA/Protein sequence:', height=200)

motif_input = st.sidebar.text_input('Enter motifs (supports comma-separated motif input [Example: ATG,TAAC,GAC], supports regex):')
mismatch_tolerance = st.sidebar.slider("Allowed mismatches:", 0, 3, 0)

# Checkbox to search in reverse complement
reverse_complement = st.sidebar.checkbox("Search in Reverse Complement")

sequence = re.sub(r'\s+', '', sequence)

motif_colors = ["#FF6347", "#FFD700", "#32CD32", "#1E90FF", "#FF69B4"]

# Compute reverse complement if needed
if reverse_complement:
    rev_complement_seq = str(Seq.Seq(sequence).reverse_complement())
else:
    rev_complement_seq = None

if st.sidebar.button('🔍 Find Motifs'):
    if not sequence or not motif_input:
        st.error('⚠️ Please enter a sequence and motifs to search for.')
    else:
        motifs = [m.strip() for m in motif_input.split(',')]
        matches_dict = {}

        # Search each motif in both the original sequence and its reverse complement (if enabled)
        for i, motif in enumerate(motifs):
            color = motif_colors[i % len(motif_colors)]
            # Exact or approximate match search in original sequence
            if mismatch_tolerance == 0:
                orig_matches = [(m.start(), m.end(), motif, "+") for m in re.finditer(motif, sequence)]
            else:
                orig_matches = [(m.start, m.end, motif, "+") for m in find_near_matches(motif, sequence, max_l_dist=mismatch_tolerance)]
            all_matches = orig_matches

            # Search in reverse complement if enabled and convert coordinates
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

        # Build highlighted sequence view using original coordinates
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
The highlighted sequence view allows you to visually inspect where the motifs occur within the sequence.  
This is useful for confirming motif locations and understanding the spatial distribution of functional elements.
""")

        # Build results table with detailed match information
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
The results table enumerates each detected motif along with its genomic coordinates and strand orientation.  
Researchers can use this tabulated data to correlate motif presence with functional genomic regions.
""")

            def generate_summary_pdf():
                pdf = FPDF()
                pdf.set_auto_page_break(auto=True, margin=15)
                pdf.add_page()
                pdf.set_font("Arial", size=12)
                pdf.cell(200, 10, "Sequence Motif Finder - Summary Report", ln=True, align='C')
                pdf.ln(10)
                pdf.cell(200, 10, f"Total Sequence Length: {len(sequence)}", ln=True)
                pdf.cell(200, 10, f"Motifs Searched: {', '.join(motifs)}", ln=True)
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
            
            st.subheader("📄 Download Summary Report of the Motif Location")
            pdf_data = generate_summary_pdf()
            st.download_button("📥 Click here to Download PDF Summary", data=BytesIO(pdf_data), file_name="Motif_Summary_Report.pdf", mime="application/pdf")
            st.markdown("""
**Scientific Explanation:**  
The downloadable PDF report provides a complete summary of your analysis, including motif parameters and match details,  
which is useful for documentation and further scientific reporting.
""")

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
This plot visualizes the spatial distribution and intensity of motif occurrences along the sequence.  
High-intensity regions may correspond to functionally important domains.
""")

            st.subheader("🔥 Motif Frequency Bar Chart")
            motif_counts = Counter([motif for motif, _ in matches_dict.items()])
            fig_bar = px.bar(x=list(motif_counts.keys()), y=list(motif_counts.values()),
                             labels={'x':'Motif', 'y':'Count'}, title='Motif Frequency')
            st.plotly_chart(fig_bar, use_container_width=True)
            st.markdown("""
**Scientific Explanation:**  
The frequency bar chart compares the occurrence of each motif.  
It helps identify which motifs are most prevalent, which may be linked to regulatory or structural functions.
""")

            st.subheader("📌 Motif Distribution Across Sequence")
            fig_line = px.line(x=list(range(len(sequence))), y=heatmap_data,
                               labels={'x':'Position', 'y':'Motif Occurrences'})
            st.plotly_chart(fig_line, use_container_width=True)
            st.markdown("""
**Scientific Explanation:**  
The line plot illustrates the overall distribution of motif occurrences along the sequence.  
It can reveal clusters or hotspots that might be critical for biological activity.
""")

            st.subheader("🧬 GC Content Visualization")
            gc_content = [(sequence[:i+1].count('G') + sequence[:i+1].count('C'))/(i+1) for i in range(len(sequence))]
            fig_gc = px.line(x=list(range(len(sequence))), y=gc_content,
                             labels={'x':'Position', 'y':'GC Content'})
            st.plotly_chart(fig_gc, use_container_width=True)
            st.markdown("""
**Scientific Explanation:**  
The GC content plot shows the proportion of guanine and cytosine nucleotides along the sequence,  
which is important for understanding the thermal stability and structural characteristics of the DNA.
""")

            st.subheader("🔥 Heatmap Visualization")
            heatmap_fig = px.imshow([heatmap_data], color_continuous_scale='viridis', aspect='auto')
            st.plotly_chart(heatmap_fig, use_container_width=True)
            st.markdown("""
**Scientific Explanation:**  
The heatmap provides a visual intensity map of motif occurrences, making it easy to identify  
regions with high or low motif activity. This is useful for pinpointing areas of potential regulatory importance.
""")

            st.subheader("🔥 Motif Co-occurrence Heatmap")
            co_matrix = pd.DataFrame([[1 if motif in sequence[start:end] else 0 for motif in motifs] 
                                      for _, start, end, _, _ in results], columns=motifs)
            fig_heat, ax = plt.subplots()
            sns.heatmap(co_matrix.corr(), annot=True, cmap='coolwarm', ax=ax)
            st.pyplot(fig_heat)
            st.markdown("""
**Scientific Explanation:**  
The co-occurrence heatmap displays the correlation between the occurrence patterns of different motifs.  
It can help identify motifs that tend to appear together, suggesting possible functional interactions.
""")

            st.subheader("🏆 Sequence Logo Plot")
            # Generate sequence logos and display them in a grid (2 per row)
            logo_figures = []
            for motif, data in matches_dict.items():
                seqs = []
                for start, end, _, strand in data['matches']:
                    if strand == "+":
                        seqs.append(sequence[start:end])
                    else:
                        seqs.append(str(Seq.Seq(sequence[start:end]).reverse_complement()))
                st.markdown(f"**Sequence Logo for motif: {motif}**")
                # Only generate a logo if at least 2 sequences are present and they have uniform length.
                if len(seqs) > 1 and len(set(len(s) for s in seqs)) == 1:
                    count_mat = logomaker.alignment_to_matrix(sequences=seqs, to_type='counts')
                    fig_logo, ax = plt.subplots(figsize=(3, 2))
                    logomaker.Logo(count_mat, ax=ax)
                    ax.set_title(f"{motif}", fontsize=10)
                    ax.set_xticks([])  # remove x-ticks for clarity
                    ax.set_yticks([])  # remove y-ticks for clarity
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
The sequence logo plot summarizes the conservation and variability of the aligned motif matches.  
- **Letter Heights:** Reflect the frequency (information content) at each position.  
- **Conserved Positions:** High letter stacks indicate critical, functionally important positions.  
For example, if variations like `ATGCGT`, `ATGAGT`, and `ATGCGT` are observed, the logo may highlight that the fourth position is variable.
""")

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
This circular (polar) visualization maps motif occurrences onto a circular representation of the sequence.  
It is especially useful for genomes or plasmids where spatial organization is critical, allowing users  
to intuitively explore motif distribution around the entire structure.
""")

            st.subheader("📊 3D Motif Distribution Plot")
            # Map each motif to a numeric index for the y-axis
            motif_index = {motif: idx for idx, motif in enumerate(motifs)}
            x_vals, y_vals, z_vals, text_vals = [], [], [], []
            for motif, data in matches_dict.items():
                for start, end, _, strand in data['matches']:
                    x_vals.append(start)  # Start position of the match
                    y_vals.append(motif_index[motif])  # Numeric index representing the motif
                    z_vals.append(end - start)  # Match length
                    text_vals.append(f"{motif} ({strand})<br>Pos: {start}-{end}")
            if x_vals:
                fig_3d = go.Figure(data=[go.Scatter3d(
                    x=x_vals,
                    y=y_vals,
                    z=z_vals,
                    mode='markers',
                    marker=dict(
                        size=5,
                        color=z_vals,  # Color scale based on match length
                        colorscale='Viridis',
                        opacity=0.8
                    ),
                    text=text_vals,
                    hoverinfo='text'
                )])
                # Update y-axis tick labels to show motif names
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
The 3D scatter plot provides an additional dimension to explore motif characteristics.  
- **X-axis:** Represents the starting position of each motif occurrence.  
- **Y-axis:** Differentiates between motif types.  
- **Z-axis:** Reflects the length of the match, which may indicate variation in motif structure.  
This multidimensional view is useful for complex datasets where spatial and quantitative relationships are important.
""")
