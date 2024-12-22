import numpy as np
from scipy.stats import chi2
import tkinter as tk
from tkinter import messagebox, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class LFSR:
    def __init__(self, polynomial, seed):
        self.polynomial = polynomial
        self.state = seed
        self.n_bits = len(polynomial)

    def step(self):
        new_bit = sum(((self.state >> i) & 1) * self.polynomial[i] for i in range(self.n_bits)) % 2
        self.state = ((self.state >> 1) | (new_bit << (self.n_bits - 1))) & ((1 << self.n_bits) - 1)
        return new_bit

    def generate(self, n):
        sequence = []
        for _ in range(n):
            sequence.append(self.step())
        return sequence

def calculate_chi_square(sequence):
    """
    Выполняет проверку качества сгенерированной последовательности с использованием критерия χ²-Пирсона.
    :param sequence: Список сгенерированных чисел.
    :return: Значение статистики χ², p-value, и результат проверки гипотезы.
    """
    N = len(sequence)
    M = max(sequence) + 1  # Максимально возможный элемент последовательности
    K = M if M < 10 else int(5 * np.log10(N))  # Количество интервалов

    # Создаем интервалы
    intervals = np.linspace(0, M, K + 1)  # Создаем K равных интервалов
    counts, _ = np.histogram(sequence, bins=intervals)  # Считаем количество элементов в каждом интервале
    v_i = counts
    P_i = 1 / K  # Теоретическая вероятность для равномерного распределения

    # Вычисляем статистику S*
    S_star = N * np.sum(((v_i / N - P_i) ** 2) / P_i)

    # Вычисляем степень свободы (r = K - 1) и p-value
    r = K - 1
    critical_value = chi2.ppf(0.95, df=r)  # Критическое значение для уровня значимости 0.05
    p_value = 1 - chi2.cdf(S_star, df=r)

    # Проверяем гипотезу
    hypothesis = "Rejected" if S_star > critical_value else "Not Rejected"

    return S_star, p_value, hypothesis, critical_value, K
def calculate_period(sequence):
    for i in range(1, len(sequence)):
        if sequence[:i] == sequence[i:i + i]:
            return i
    return len(sequence)


def chi_square_test(sequence, num_bins=2):
    unique, counts = np.unique(sequence, return_counts=True)
    observed = counts
    expected = [len(sequence) / num_bins] * num_bins
    return chi2(observed, expected)


def run_simulation():
    try:
        # Получение параметров из интерфейса
        polynomial_str = polynomial_var.get()
        seed_str = seed_entry.get()
        n = int(n_entry.get())

        polynomial = polynomial_options.get(polynomial_str)
        if not polynomial:
            raise ValueError("Invalid polynomial selected.")
        if not seed_str.startswith("0b"):
            raise ValueError("Seed must be in binary format (e.g., 0b1010101).")

        seed = int(seed_str, 2)
        if seed <= 0 or seed >= (1 << len(polynomial)):
            raise ValueError(f"Seed must be a binary number with {len(polynomial)} bits.")

        # Генерация последовательности
        lfsr = LFSR(polynomial, seed)
        sequence = lfsr.generate(n)
        period = calculate_period(sequence)
        S_star, p_value, hypothesis, critical_value, K = calculate_chi_square(sequence)

        # Обновление текстового поля
        result_text.delete(1.0, tk.END)
        result_text.insert(tk.END, f'Polynomial: {polynomial_str}\n')
        result_text.insert(tk.END, f'Seed: {seed_str}\n')
        result_text.insert(tk.END, f'Number of bits (n): {n}\n')
        result_text.insert(tk.END, f'Sequence: {"".join(map(str, sequence))}\n')
        result_text.insert(tk.END, f'Period: {period}\n')
        result_text.insert(tk.END, f'Chi-square statistic (S*): {S_star:.4f}\n')
        result_text.insert(tk.END, f'p-value: {p_value:.4f}\n')
        result_text.insert(tk.END, f'Critical value: {critical_value:.4f}\n')
        result_text.insert(tk.END, f'Number of intervals (K): {K}\n')
        result_text.insert(tk.END, f'Hypothesis test result: {hypothesis}\n')

        # Отображение графика
        display_graph(sequence, period, S_star, p_value, hypothesis)

    except Exception as e:
        messagebox.showerror("Error", str(e))


def display_graph(sequence, period, S_star, p_value, hypothesis):
    # Построение графика последовательности
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(sequence, drawstyle="steps-pre", marker="o", label="Sequence")
    ax.set_title(f"LFSR Sequence Visualization\nPeriod: {period} | S*: {S_star:.2f}, p-value: {p_value:.4f}\nHypothesis: {hypothesis}")
    ax.set_xlabel("Step")
    ax.set_ylabel("Bit Value")
    ax.legend()

    # Вставка графика в Tkinter
    for widget in graph_frame.winfo_children():
        widget.destroy()
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()


# Создание GUI
root = tk.Tk()
root.title("LFSR Generator with Visualization")

# Параметры многочлена
polynomial_options = {
    "f(x) = x^7 + x + 1": [1, 0, 0, 0, 0, 0, 1, 1],
    "f(x) = x^7 + x^5 + x^3 + 1": [1, 0, 1, 0, 1, 0, 0, 1],
    "f(x) = x^7 + x^6 + x^5 + x^2 + 1": [1, 1, 1, 0, 0, 1, 0, 1]
}

# Выбор многочлена
polynomial_var = tk.StringVar(value=list(polynomial_options.keys())[0])
tk.Label(root, text="Choose Polynomial:").pack(pady=5)
polynomial_menu = tk.OptionMenu(root, polynomial_var, *polynomial_options.keys())
polynomial_menu.pack(pady=5)

# Поле для ввода начального состояния
tk.Label(root, text="Seed (e.g., 0b1010101):").pack(pady=5)
seed_entry = tk.Entry(root)
seed_entry.pack(pady=5)

# Поле для ввода длины последовательности
tk.Label(root, text="Number of bits (n):").pack(pady=5)
n_entry = tk.Entry(root)
n_entry.pack(pady=5)

# Кнопка запуска симуляции
run_button = tk.Button(root, text="Run Simulation", command=run_simulation)
run_button.pack(pady=20)

# Текстовое поле для отображения результатов
result_text = scrolledtext.ScrolledText(root, width=50, height=15)
result_text.pack(pady=10)

# Фрейм для графика
graph_frame = tk.Frame(root)
graph_frame.pack(pady=10)

root.mainloop()
