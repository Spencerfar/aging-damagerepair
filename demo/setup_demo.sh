git clone https://github.com/SinclairLab/frailty
python convert_excel.py
mkdir -p ../datasets/
python ../clean_data/clean_schultz.py --folder './'
