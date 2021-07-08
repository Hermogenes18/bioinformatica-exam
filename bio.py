from flask import Flask, render_template, request
import aplicaciones

app = Flask(__name__)

@app.route('/')
def student():
   return render_template('index.html')

@app.route('/metodos',methods = ['POST', 'GET'])
def login():
   return render_template()

if __name__ == '__main__':
   app.run(debug = True)