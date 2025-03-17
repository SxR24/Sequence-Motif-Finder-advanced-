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
from Bio import SeqIO
import io

# Streamlit Page Configuration
st.set_page_config(page_title="Sequence Motif Finder", layout="wide")
st.title('🔬 Advanced Sequence Motif Finder')

# Sidebar Inputs
st.sidebar.header('⚙️ Input Options')

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

motif_input = st.sidebar.text_input('Enter motifs (comma-separated, supports regex):')
mismatch_tolerance = st.sidebar.slider("Allowed mismatches:", 0, 3, 0)

sequence = re.sub(r'\s+', '', sequence)

motif_colors = ["#FF6347", "#FFD700", "#32CD32", "#1E90FF", "#FF69B4"]

if st.sidebar.button('🔍 Find Motifs'):
    if not sequence or not motif_input:
        st.error('⚠️ Please enter a sequence and motifs to search for.')
    else:
        motifs = [m.strip() for m in motif_input.split(',')]
        matches_dict = {}

        for i, motif in enumerate(motifs):
            color = motif_colors[i % len(motif_colors)]
            if mismatch_tolerance == 0:
                matches = [(m.start(), m.end(), motif) for m in re.finditer(motif, sequence)]
            else:
                matches = [(m.start, m.end, motif) for m in find_near_matches(motif, sequence, max_l_dist=mismatch_tolerance)]
            
            if matches:
                matches_dict[motif] = {'matches': matches, 'color': color}

        positions = []
        for motif, data in matches_dict.items():
            color = data['color']
            for start, end, _ in data['matches']:
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

        results = []
        for motif, data in matches_dict.items():
            for start, end, _ in data['matches']:
                results.append([motif, start, end, sequence[start:end]])
        
        if results:
            df_results = pd.DataFrame(results, columns=['Motif', 'Start', 'End', 'Matched Sequence'])
            st.subheader('🔍 Motif Matches')
            st.dataframe(df_results)

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
                    for start, end, _ in data['matches']:
                        pdf.cell(200, 10, f"Motif: {motif}, Start: {start}, End: {end}, Sequence: {sequence[start:end]}", ln=True)
                return pdf.output(dest='S').encode('latin1')
            
            st.subheader("📄 Download Summary Report of the motif Location")
            pdf_data = generate_summary_pdf()
            st.download_button("📥 Click here to Download PDF Summary", data=BytesIO(pdf_data), file_name="Motif_Summary_Report.pdf", mime="application/pdf")

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
                for start, end, _ in data['matches']:
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

            
            
            st.subheader("🔥 Motif Frequency Bar Chart")
            motif_counts = Counter([motif for motif, _ in matches_dict.items()])
            fig_bar = px.bar(x=motif_counts.keys(), y=motif_counts.values(), labels={'x':'Motif', 'y':'Count'}, title='Motif Frequency')
            st.plotly_chart(fig_bar, use_container_width=True)

            st.subheader("📌 Motif Distribution Across Sequence")
            fig_line = px.line(x=list(range(len(sequence))), y=heatmap_data, labels={'x':'Position', 'y':'Motif Occurrences'})
            st.plotly_chart(fig_line, use_container_width=True)

            st.subheader("🧬 GC Content Visualization")
            gc_content = [(sequence[:i+1].count('G') + sequence[:i+1].count('C')) / (i+1) for i in range(len(sequence))]
            fig_gc = px.line(x=list(range(len(sequence))), y=gc_content, labels={'x':'Position', 'y':'GC Content'})
            st.plotly_chart(fig_gc, use_container_width=True)

            # Heatmap visualization
            st.subheader("🔥 Heatmap Visualization")
            heatmap_fig = px.imshow([heatmap_data], color_continuous_scale='viridis', aspect='auto')
            st.plotly_chart(heatmap_fig, use_container_width=True)

            st.subheader("🔥 Motif Co-occurrence Heatmap")
            co_matrix = pd.DataFrame([[1 if motif in sequence[start:end] else 0 for motif in motifs] for _, start, end, _ in results], columns=motifs)
            fig_heat = sns.heatmap(co_matrix.corr(), annot=True, cmap='coolwarm')
            fig, ax = plt.subplots()
            sns.heatmap(co_matrix.corr(), annot=True, cmap='coolwarm', ax=ax)
            st.pyplot(fig)

            st.subheader("🏆 Sequence Logo Plot (Placeholder)")
            st.write("Sequence Logo visualization will be added soon.")

            st.subheader("🌀 Circular Genome Visualization")
            features = [GraphicFeature(start=start, end=end, strand=+1, color=motif_colors[i % len(motif_colors)], label=motif) for i, (motif, data) in enumerate(matches_dict.items()) for start, end, _ in data['matches']]
            record = GraphicRecord(sequence_length=len(sequence), features=features)
            fig, ax = plt.subplots(figsize=(10, 2))
            record.plot(ax=ax)
            st.pyplot(fig)
