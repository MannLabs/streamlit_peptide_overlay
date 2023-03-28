import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import math
import os 

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def get_overlap_positions(peptide_positions):
    overlap_positions = []
    for i, (_, pos1, *_) in enumerate(peptide_positions[:-1]):
        for j, (_, pos2, *_) in enumerate(peptide_positions[i + 1:], i + 1):
            if pos1 <= pos2 < pos1 + len(peptide_positions[i][0]) or pos2 <= pos1 < pos2 + len(peptide_positions[j][0]):
                overlap_positions.append((i, j))
    return overlap_positions

def find_all_positions(string, sub_string):
    start = 0
    positions = []
    while True:
        index = string.find(sub_string, start)
        if index == -1:
            break
        positions.append((sub_string, index))
        start = index + 1
    return positions

st.set_option("deprecation.showfileUploaderEncoding", False)

st.title("Peptide Visualization")

protein_sequence = st.text_area("Enter the protein sequence:")
peptide_input = st.text_area("Enter the digested peptides separated by commas:")

c1, c2 = st.columns(2)
offset = c1.number_input('Offset for overlapping peptides', min_value=0.001, max_value=1.0, value=0.1, step=0.01)
offset_file = c2.number_input('Offset for groups', min_value=0.001, max_value=1.0, value=0.5, step=0.01)

uploaded_files = st.file_uploader("Upload one or several CSV file (Should contain a Header with Sequence and Score column)", type=["csv"], accept_multiple_files=True)

if protein_sequence and peptide_input and uploaded_files is not None:

    peptides = peptide_input.split(",")

    peptide_positions = []
    for peptide in peptides:

        pos = find_all_positions(protein_sequence, peptide)

        if len(pos) == 0:
            st.warning(f"In-silico peptide {peptide} not found in protein sequence")
        else:
            peptide_positions.extend(pos)

    overlap_positions = get_overlap_positions(peptide_positions)

    fig = go.Figure()

    y_offset = {}
    for i, (peptide, position) in enumerate(peptide_positions):
        y_value = 0
        for a, b in overlap_positions:
            if i == b:
                y_value += offset

        y_offset[i] = y_value

    for i, (peptide, position) in enumerate(peptide_positions):
        fig.add_trace(
            go.Scatter(
                x=list(range(position, position + len(peptide))),
                y=[y_offset[i]] * len(peptide),
                mode="lines",
                marker=dict(size=10, color='black'),
                name="Digested Peptides",
                legendgroup="Input Peptides",
                hovertemplate=(
                f"Peptide: {peptide}<br>"
                f"Start: {position}<br>"
                f"End: {position+len(peptide)}<br>"
                f"Length: {len(peptide)}<br>"),
                showlegend=i == 0,
            )
        )


    off = math.ceil(max(y_offset.values()))

    for file_idx, uploaded_file in enumerate(uploaded_files):
        df = pd.read_csv(uploaded_file)

        st.write(f"Uploaded Data {uploaded_file.name}:")
        st.write(df)

        df_peptide_positions = []
        for idx, row in df.iterrows():
            peptide = row["Sequence"]
            pos = find_all_positions(protein_sequence, peptide)

            if len(pos) == 0:
                st.warning(f"Peptide {peptide} not found in protein sequence")
            else:
                df_peptide_positions.extend(pos)

        df_peptide_positions.sort(key=lambda x: x[2], reverse=True)
        df_overlap_positions = get_overlap_positions(df_peptide_positions)

        df_y_offset = {}
        for i, (peptide, position, _) in enumerate(df_peptide_positions):
            y_value = off+offset_file
            for a, b in df_overlap_positions:
                if i == b:
                    y_value += offset

            df_y_offset[i] = y_value
        
        off = math.ceil(max(df_y_offset.values()))
        fname = uploaded_file.name
        fname_short = os.path.splitext(fname)[0]
        for i, (peptide, position, _) in enumerate(df_peptide_positions):
            
            score = df_peptide_positions[i][2]
            fig.add_trace(
                go.Scatter(
                    x=list(range(position, position + len(peptide))),
                    y=[df_y_offset[i]] * len(peptide),
                    mode="lines",
                    marker=dict(size=10, color=colors[file_idx]),
                    name=f'{fname_short}', # name=peptide,
                    hovertemplate=(
                    f"Peptide: {peptide}<br>"
                    f"Start: {position}<br>"
                    f"End: {position+len(peptide)}<br>"
                    f"Length: {len(peptide)}<br>"
                    f"Score: {score}<br>"),
                    legendgroup="CSV Peptides",
                    showlegend=i == 0,
                )
            )


    fig.update_xaxes(title="Protein Sequence Position")
    #fig.update_yaxes(title="Peptides", range=[0, max_y_value + 1], visible=False, showticklabels=False)
    fig.update_layout(title="Peptide Positions in the Protein Sequence")

    st.plotly_chart(fig)