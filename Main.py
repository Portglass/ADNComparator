import ADN_library as al
import Visual_Streamlit as vs
import streamlit as st
import matplotlib

matplotlib.use("Agg")
from Bio import SeqIO
from io import StringIO
import altair as alt


def main():
    choice = vs.head()
    list_seq = []

    if choice == 'Intro':
        st.subheader("Intro")

    elif choice == 'DNA':
        st.subheader("DNA details")
        seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"], accept_multiple_files=True)
        for file in seq_file:
            st.subheader(file.name)
            stringio = StringIO(file.getvalue().decode("utf-8"))
            sequence = ""
            for record in SeqIO.parse(stringio, 'fasta'):
                sequence = str(record.seq)
                st.write(f'Length of sequence: {len(sequence)}')
            list_seq.append(list(sequence))
            df_ADN = al.fromStringToDataframe(list(sequence))
            vs.bar_chart(df_ADN)
            _, df_pourcent_ADN = al.pourcentNucleo(list(sequence))
            bar_chart = alt.Chart(df_pourcent_ADN).mark_bar().encode(
                y='Pourcent:Q',
                x='Nucleotide:O',
            )
            st.altair_chart(bar_chart, use_container_width=True)



    elif choice == "Comparaison entre deux ADN":#Fonctionne pas
        st.subheader("Comparaison entre deux ADN")
        print(len(list_seq))
        if len(list_seq) <= 1:
            st.subheader("Aucun ADN n'a été enregistrer veuillez vous rendre dans l'onglet DNA")
        else:
            st.write("Le nombre de mistmatch est " + str(al.nbMismatch(list_seq[0], list_seq[1])))


    elif choice == 'About':
        st.subheader("About")
        st.write("Réalisé par Martin R2villon, étudiant à l'EBI et Efrei")


if __name__ == '__main__':
    main()
