import tkinter as tk
#from tkinter.ttk import *
from tkinter import ttk

from Powell_Gauss_Seidel import powell, create_ksi, get_value_from_function
from GaussSeidel import gauss_siedel
import sympy
import matplotlib.pyplot as plt
import numpy as np
import tkinter.font as font
#------------------------------------------------------
root = tk.Tk()
root.wm_title(' Powell - Gauss - Seidel ')
root.geometry('690x475')
#------------------------------------------------------
s = ttk.Style()
s.configure('Red.TLabelframe.Label', font=('courier', 15, 'bold'))
s.configure('Red.TLabelframe.Label', foreground='blue')
labelframe = ttk.LabelFrame(root, text = " 1. Funkcja celu: ", style="Red.TLabelframe")
labelframe.grid(row=0, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)

left = tk.Label(labelframe, text="Wybierz funkcje:", font='Verdana 10 bold')
left.grid(row=4, column=0, sticky='WE', padx=5, pady=5, ipadx=5, ipady=5)

combo = ttk.Combobox(labelframe, width=40, font="Verdana 10")
combo['values'] = ('', '(x1-2)**2 + (x2-1)**2', 'x1**4 + x2**4 - x1**2 - x2**2', '(x1 - 2)**2 + (x1-x2**2)**2')#, 'x1**4 + x2**4 - 2*x1**2*x2 - 4*x1 + 3')
combo.current(0)
combo.grid(row=5, column=0, sticky='WE', padx=5, pady=5, ipadx=5, ipady=5)
#------------------------------------------------------
labelframe2 = ttk.Labelframe(root, text=" 2. Podaj punkt początkowy: ", style="Red.TLabelframe")
labelframe2.grid(row=0, column=8, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)

left1 = tk.Label(labelframe2, text="Wartości:", font='Verdana 10 bold')
left1.grid(row=1, column=0, sticky='WE', padx=5, pady=5, ipadx=5, ipady=5)

lbl_x1 = tk.Label(labelframe2, text="X : ", font='Verdana 10 bold')
lbl_x1.grid(row=2, column=0, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)
txt_x1 = tk.Entry(labelframe2, width=30, font='Verdana 10 bold')
txt_x1.grid(row=2, column=1, columnspan=5, sticky='E', padx=5, pady=5, ipadx=5, ipady=5)
#------------------------------------------------------
labelframe3 = ttk.Labelframe(root, text=" 3. Ograniczenia: ", style="Red.TLabelframe")
labelframe3.grid(row=3, column=0, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)

txt_g1 = tk.Text(labelframe3, height=7, width=37)
txt_g1.grid(row=4, column=2, sticky='WE', pady=15)
lbl_g1 = tk.Label(labelframe3, text="<= 0", font='Helvetica 18 bold')
lbl_g1.grid(row=4, column=10, sticky='WE', pady=15)
#------------------------------------------------------
labelframe4 = ttk.Labelframe(root, text=" 4. Kryteria stopu: ", style="Red.TLabelframe")
labelframe4.grid(row=3, column=8, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)

epsilon1 = tk.Label(labelframe4, text="epsilon1 : ", font='Verdana 10 bold')
epsilon1.grid(row=4, column=10, sticky='WE', padx=5, pady=3)
txt_epsilon1 = tk.Entry(labelframe4, width=31, font='Verdana 10 bold')
txt_epsilon1.insert(tk.END, '0.1')
txt_epsilon1.grid(row=4, column=17, sticky='WE', padx=5, pady=3)

epsilon2 = tk.Label(labelframe4, text="epsilon2 : ", font='Verdana 10 bold')
epsilon2.grid(row=5, column=10, sticky='WE', padx=5, pady=3)
txt_epsilon2 = tk.Entry(labelframe4, width=18, font='Verdana 10 bold')
txt_epsilon2.insert(tk.END, '0.001')
txt_epsilon2.grid(row=5, column=17, sticky='WE', padx=5, pady=3)

epsilon3 = tk.Label(labelframe4, text="epsilon3 : ", font='Verdana 10 bold')
epsilon3.grid(row=6, column=10, sticky='WE', padx=5, pady=3)
txt_epsilon3 = tk.Entry(labelframe4, width=18, font='Verdana 10 bold')
txt_epsilon3.insert(tk.END, '0.001')
txt_epsilon3.grid(row=6, column=17, sticky='WE', padx=5, pady=3)

epsilon4 = tk.Label(labelframe4, text="E (krok) : ", font='Verdana 10 bold')
epsilon4.grid(row=7, column=10, sticky='WE', padx=5, pady=3)
txt_epsilon4 = tk.Entry(labelframe4, width=18, font='Verdana 10 bold')
txt_epsilon4.insert(tk.END, '5')
txt_epsilon4.grid(row=7, column=17, sticky='WE', padx=5, pady=3)

epsilon5 = tk.Label(labelframe4, text="L : ", font='Verdana 10 bold')
epsilon5.grid(row=8, column=10, sticky='WE', padx=5, pady=3)
txt_epsilon5 = tk.Entry(labelframe4, width=18, font='Verdana 10 bold')
txt_epsilon5.insert(tk.END, '1000')
txt_epsilon5.grid(row=8, column=17, sticky='WE', padx=5, pady=3)
#------------------------------------------------------
labelframe5 = ttk.Labelframe(root, text=" 5. Kroki algorytmu: ", style="Red.TLabelframe")
labelframe5.grid(row=8, column=0, sticky='W', padx=5, pady=5, ipadx=5, ipady=5)

txt_powell = tk.Text(labelframe5, height=10, width=43, font='Verdana 10 bold')
txt_powell.grid(row=9, column=0, sticky='WE', padx=3, pady=3)
#------------------------------------------------------


def calc():
    x1, x2, x3, x4, x5 = sympy.symbols("x1 x2 x3 x4 x5")
    f = sympy.sympify(combo.get())

    g = []
    print(txt_x1.get())
    for line in txt_g1.get('1.0', 'end-1c').splitlines():
        if line:
            g.append(sympy.sympify(line))

    c_min = float(txt_epsilon1.get())
    epsilon_j = min(float(txt_epsilon2.get()), float(txt_epsilon3.get()))
    e_0 = int(txt_epsilon4.get())
    L = int(txt_epsilon5.get())

    x_str = list(txt_x1.get().split(';'))
    x = list(map(int, x_str))

    ksi = create_ksi(len(x))
    c = 5

    steps_x1 = []
    steps_x2 = []

    if not g:
        steps_x1, steps_x2 = gauss_siedel(f, x, ksi, e_0, epsilon_j, L)

        txt_powell.delete('1.0', tk.END)
        txt_powell.insert(tk.END, 'Gauss-Seidel')
        txt_powell.insert(tk.END, '\nx = ({0},{1})\nF(x) = {2}'.format(steps_x1[0], steps_x2[0], get_value_from_function(f, [steps_x1[0], steps_x2[0]])))

        for i in range(len(steps_x1) - 1):
            txt_powell.insert(tk.END, '\nx = ({0},{1})\nF(x) = {2}'.format(steps_x1[i+1], steps_x2[i+1], get_value_from_function(f, [steps_x1[i+1], steps_x2[i+1]])))
    else:
        powell_steps, c_steps = powell(f, g, x, c, c_min, e_0, epsilon_j, L, ksi)
        print(powell_steps, c_steps)

        txt_powell.delete('1.0', tk.END)
        txt_powell.insert(tk.END, 'Powell')
        txt_powell.insert(tk.END, '\nx = ({0},{1})\nF(x) = {2}, c = {3}'.format(powell_steps[0][0],
                                                                      powell_steps[0][1],
                                                                      get_value_from_function(f, [powell_steps[0][0], powell_steps[0][1]]), c_steps[0]))

        for i in range(len(powell_steps) - 1):
            txt_powell.insert(tk.END, '\nx = ({0},{1})\nF(x) = {2}, c= {3}'.format(powell_steps[i+1][0],
                                                                         powell_steps[i+1][1],
                                                                         get_value_from_function(f, [powell_steps[i+1][0], powell_steps[i+1][1]]), c_steps[i+1]))


def draw():
    x1, x2, x3, x4, x5 = sympy.symbols("x1 x2 x3 x4 x5")
    f = sympy.sympify(combo.get())

    g = []
    print(txt_x1.get())
    for line in txt_g1.get('1.0', 'end-1c').splitlines():
        if line:
            g.append(sympy.sympify(line))

    c_min = float(txt_epsilon1.get())
    epsilon_j = min(float(txt_epsilon2.get()), float(txt_epsilon3.get()))
    e_0 = int(txt_epsilon4.get())
    L = int(txt_epsilon5.get())

    x_str = list(txt_x1.get().split(';'))
    x = list(map(int, x_str))

    ksi = create_ksi(len(x))
    c = 5

    if len(x) == 2:
        if not g:
            steps_x1, steps_x2 = gauss_siedel(f, x, ksi, e_0, epsilon_j, L)

            x_begin = min(steps_x1[0], steps_x1[-1])
            x_end = max(steps_x1[0], steps_x1[-1])

            x_begin = x_begin - 3
            x_end = x_end + 3

            y_begin = min(steps_x2[0], steps_x2[-1])
            y_end = max(steps_x2[0], steps_x2[-1])

            y_begin = y_begin - 3
            y_end = y_end + 3

            x = np.linspace(int(x_begin), int(x_end), 50)
            y = np.linspace(int(y_begin), int(y_end), 50)

            X, Y = np.meshgrid(x, y)
            Z = np.zeros((50, 50))

            for i in range(50):
                for j in range(50):
                    Z[i][j] = get_value_from_function(f, [X[i][j], Y[i][j]])

            plt.contourf(X, Y, Z, 100, cmap='RdGy')
            plt.colorbar()
            for i in range(len(steps_x1) - 1):
                plt.plot([steps_x1[i], steps_x1[i + 1]], [steps_x2[i], steps_x2[i + 1]], marker=".", markersize=10, color="green")

            plt.show()
        else:
            powell_steps, c_steps = powell(f, g, x, c, c_min, e_0, epsilon_j, L, ksi)

            x_begin = min(powell_steps[0][0], powell_steps[-1][0])
            x_end = max(powell_steps[0][0], powell_steps[-1][0])

            x_begin = x_begin - 3
            x_end = x_end + 3

            y_begin = min(powell_steps[0][1], powell_steps[-1][1])
            y_end = max(powell_steps[0][1], powell_steps[-1][1])

            y_begin = y_begin - 3
            y_end = y_end + 3

            x = np.linspace(int(x_begin), int(x_end), 50)
            y = np.linspace(int(y_begin), int(y_end), 50)

            X, Y = np.meshgrid(x, y)
            print(X)
            Z = np.zeros((50, 50))

            for i in range(50):
                for j in range(50):
                    Z[i][j] = get_value_from_function(f, [X[i][j], Y[i][j]])

            plt.contourf(X, Y, Z, 100, cmap='RdGy')
            plt.colorbar()
            for i in range(len(powell_steps) - 1):
                plt.plot([powell_steps[i][0], powell_steps[i + 1][0]], [powell_steps[i][1], powell_steps[i + 1][1]], marker=".", markersize=10, color="green")

            g_colors = ["black", "blue", "purple", "yellow"]
            for gi in range(len(g)):
                g_contour = np.zeros((50, 50))
                for i in range(50):
                    for j in range(50):
                        g_contour[i][j] = get_value_from_function(g[gi], [X[i][j], Y[i][j]])
                plt.contour(X, Y, g_contour, 0, colors=g_colors[gi])

            plt.show()


def clear_text():
    txt_powell.delete('1.0', tk.END)


labelframe6 = ttk.Labelframe(root, text=" 6. Rysuj: ", style="Red.TLabelframe")
labelframe6.grid(row=8, column=8, sticky='NEW')

myFont = font.Font(family='Verdana', size=10, weight='bold')

btn = tk.Button(labelframe6, height=2, width=13, text="Policz", command=calc, bg='white', fg='black', padx=7, pady=6)
btn['font'] = myFont
btn.grid(row=1, column=0)

btn2 = tk.Button(labelframe6, height=2, width=13, text="Rysuj", command=draw, bg='white', fg='black', padx=7, pady=6)
btn2['font'] = myFont
btn2.grid(row=1, column=1)

btn3 = tk.Button(labelframe6, height=2, width=13, text="Wyczyść", command=clear_text, bg='white', fg='black', padx=7, pady=6)
btn3['font'] = myFont
btn3.grid(row=1, column=2)
#------------------------------------------------------
root.mainloop()
#------------------------------------------------------
