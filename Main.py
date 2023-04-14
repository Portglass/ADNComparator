import ADN_library as al
import streamlit as st
import Visual_Streamlit as vs
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np
from PIL import Image



def main():
    choice = vs.head()

    if choice == 'Intro':
        st.subheader("Intro")

    elif choice == 'DNA':
        st.subheader("DNA details")
        seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
        if seq_file is not None:
                dna_record = SeqIO.read(seq_file, "fasta")
                # st.write(dna_record)
                dna_seq = dna_record.seq

                details = st.radio("Details", ("Description", "Sequence"))
                if details == "Description":
                    st.write(dna_record.description)
                elif details == "Sequence":
                    st.write(dna_record.seq)

    elif choice == 'About':
        st.subheader("About")


        """
        header, ADN_list = al.read_file(uploaded_file)
        df_ADN = al.fromFileToDataframe(uploaded_file)

        vs.bar_chart(df_ADN)
        _, df_pourcent_ADN = al.pourcentNucleo(ADN_list)
        bar_chart = alt.Chart(df_pourcent_ADN).mark_bar().encode(
            y='Pourcent:Q',
            x='Nucleotide:O',
        )

        st.altair_chart(bar_chart, use_container_width=True)
        """

if __name__ == '__main__':
	main()



