from flask import Flask, render_template, request
import pandas as pd
from ali_local import funcion
import blast as blst
import ali_global as ag
import jukes_cantor as jukes
<<<<<<< HEAD
from kimura import K2Pdistance
=======
import kimura as kim
import Tamura as tam
import Tajima as taj
>>>>>>> 5b151472c64980f1854c7e54c64771070b0a57fc

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import SeqIO
import os
from flask import flash, redirect, url_for
from werkzeug.utils import secure_filename

APP_ROOT = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)

@app.route('/')
def main():
    return render_template('app.html')


def read_fasta(filename):
    sequences = SeqIO.parse(filename, "fasta")
    for record in sequences:
        data = str(record.seq.upper())
    return data[0:len(data)//4]

@app.route('/send', methods=['POST'])
def send(sum=sum):
    if request.method == 'POST':
        if ('file' in request.files == True ) or ('file1' in request.files == True):
            print ("archivo")
            file = request.files['file']
            target = os.path.join(APP_ROOT, 'uploads/')
            filename = secure_filename(file.filename)
            destination = "/".join([target, filename])
            file.save(destination)
            num1 = read_fasta(destination)
            if 'file1' in request.files:
                file1 = request.files['file1']
                filename1 = secure_filename(file1.filename)
                destination1 = "/".join([target, filename1])
                file1.save(destination1)
                num2 = read_fasta(destination1)
            else:
                num2 = request.form['num2']

        else:
            num1 = request.form['num1']
            num2 = request.form['num2']
            print("IF")
        operation = request.form['operation']

        if operation == 'add':
            #sum = float(num1) + float(num2)
  
            values,A,B = funcion(str(num1),str(num2),d=-5,m_b=False)
            return render_template('app.html', sum=values[0], seq=values[1], total=values[2])

        elif operation == 'subtract':
            #sum = float(num1) - float(num2)
            values,A,B = ag.funcion(str(num1),str(num2),d=-5,m_b=False)
            return render_template('app.html', sum=values[0], seq=values[1], total=values[2])

        elif operation == 'multiply':
            value = blst.main(str(num1))
            values = value[0]
            print(values)
            print(type(value),value)
            return render_template('app.html', sum=values[0], seq=values[1], total=values[2])

        elif operation == 'divide':
            #sum = float(num1) / float(num2)
            muscle_exe = r"G:/bioinformatica/SEGUNDOEXAM ESTESI/bioinformatica-exam/muscle3.8.31_i86win32.exe"
            #print(muscle_exe, type(muscle_exe))
            in_file = destination
            out_file = os.path.join(APP_ROOT, 'aligned.fasta')
            cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
            stdout, stderr = cline(in_file)
            align = AlignIO.read(out_file, "fasta")
            for record in align:
                data = [record.id,record.seq]
            return render_template('app.html', sum=data[0],seq=data[1])
        elif operation == 'jukes':
            matrix = jukes.main(destination)
            for i in matrix:
                for j in i:
                    print(j)
            return render_template('app.html', sum=0, data=matrix)
<<<<<<< HEAD

        elif operation == 'kimura':
            values = K2Pdistance(str(num1),str(num2))
            return render_template('app.html', sum=values)
=======
        elif operation == 'kimura':
            print ("kimura")
            respuesta = kim.K2Pdistance(str(num1),str(num2))
            return render_template('app.html',total=respuesta)
        elif operation == 'tamura':
            print ("tamura")
            respuesta = tam.Tamuradistance(str(num1),str(num2))
            return render_template('app.html',total=respuesta)
        elif operation == 'tajima':
            print ("tajima")
            respuesta = taj.TNdistance(str(num1),str(num2))
            return render_template('app.html',total=respuesta)



>>>>>>> 5b151472c64980f1854c7e54c64771070b0a57fc
        else:
            return render_template('app.html')


if __name__ == "__main__":
    app.run(port=5000, debug=True)
