# bioinformatics-seq-analysis
Analysis of an unknown sequence.

## Instalar

Crear el entorno virtual con `virtualenv`:
```
virtualenv .env
o
python3 -m venv "virtualenv"
```

Activar el entorno:
```
source .env/bin/activate
```

Instalar dependencias:
```
pip install -r requirements.txt
```

### Ejemplos de Llamados

Calculo de ORF:
```
python main.py -e 1 -gb sequences/gb_unknown_seq.gb
```
Analisis del motif:
```
cd emboss
patmatmotifs -sequence orf2_protein.fasta -full
```