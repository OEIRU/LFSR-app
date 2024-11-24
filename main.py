import numpy as np
from scipy.stats import chisquare
import tkinter as tk
from tkinter import messagebox, scrolledtext


class LFSR:
    def __init__(self, polynomial, seed):
        self.polynomial = polynomial
        self.state = seed

    def step(self):
        new_bit = sum((self.state >> i) & 1 for i in range(len(self.polynomial))) % 2
        self.state = (self.state >> 1) | (new_bit << (len(self.polynomial) - 1))
        return new_bit

    def generate(self, n):
        sequence = []
        for _ in range(n):
            sequence.append(self.step())
        return sequence


def calculate_period(sequence):
    seen = {}
    for i, num in enumerate(sequence):
        if num in seen:
            return i - seen[num]
        seen[num] = i
    return len(sequence)


def chi_square_test(sequence, num_bins=2):
    observed, _ = np.histogram(sequence, bins=num_bins)
    expected = [len(sequence) / num_bins] * num_bins
    return chisquare(observed, expected)


def save_results(sequence, period, chi2_stat, p_value):
    with open('output.txt', 'w') as f:
        f.write(f'Sequence: {sequence}\n')
        f.write(f'Period: {period}\n')
        f.write(f'Chi-square statistic: {chi2_stat}, p-value: {p_value}\n')


def run_simulation():
    try:
        polynomial_str = polynomial_var.get()
        seed_str = seed_entry.get()
        n = int(n_entry.get())

        # Преобразование строки многочлена в список коэффициентов
        polynomial = polynomial_options[polynomial_str]

        # Проверка корректности ввода начального состояния
        if not seed_str.startswith("0b"):
            raise ValueError("Seed must be in binary format (e.g., 0b1010101).")

        seed = int(seed_str, 2)

        lfsr = LFSR(polynomial, seed)
        sequence = lfsr.generate(n)

        period = calculate_period(sequence)
        chi2_stat, p_value = chi_square_test(sequence)

        save_results(sequence, period, chi2_stat, p_value)

        # Отображение результатов в текстовом поле
        result_text.delete(1.0, tk.END)
        result_text.insert(tk.END, f'Sequence: {sequence}\n')
        result_text.insert(tk.END, f'Period: {period}\n')
        result_text.insert(tk.END, f'Chi-square statistic: {chi2_stat}, p-value: {p_value}\n')

    except Exception as e:
        messagebox.showerror("Error", str(e))


# GUI
root = tk.Tk()
root.title("LFSR Random Number Generator")

# Выбор многочлена
polynomial_options = {
    "f(x) = x^7 + x^3 + x^2 + x + 1": [1, 0, 0, 1, 1, 0, 1],
    "f(x) = x^7 + x^5 + x^3 + 1": [1, 0, 1, 0, 1, 0, 1],
    "f(x) = x^7 + x^6 + x^5 + x^2 + 1": [1, 1, 1, 0, 0, 1, 1]
}

polynomial_var = tk.StringVar(value=list(polynomial_options.keys())[0])
tk.Label(root, text="Choose Polynomial:").pack(pady=5)
polynomial_menu = tk.OptionMenu(root, polynomial_var, *polynomial_options.keys())
polynomial_menu.pack(pady=5)

# Поле для ввода начального состояния
tk.Label(root, text="Seed (e.g., 0b1010101):").pack(pady=5)
seed_entry = tk.Entry(root)
seed_entry.pack(pady=5)

# Поле для ввода количества бит
tk.Label(root, text="Number of bits (n):").pack(pady=5)
n_entry = tk.Entry(root)
n_entry.pack(pady=5)

# Кнопка запуска симуляции
run_button = tk.Button(root, text="Run Simulation", command=run_simulation)
run_button.pack(pady=20)

# Текстовое поле для отображения результатов
result_text = scrolledtext.ScrolledText(root, width=50, height=15)
result_text.pack(pady=10)

root.mainloop()
