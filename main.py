import numpy as np
from scipy.stats import chi2
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt


class LFSRApp:
    def __init__(self, root):
        self.root = root
        self.root.title("LFSR Generator with Visualization")

        # Interface setup
        self.setup_ui()

    def setup_ui(self):
        frame = ttk.Frame(self.root, padding="10")
        frame.grid(row=0, column=0, sticky="NSEW")

        # Input file selection
        ttk.Label(frame, text="Input File (input.txt):").grid(row=0, column=0, sticky="W")
        ttk.Button(frame, text="Run Simulation from File", command=self.run_simulation_from_file).grid(row=0, column=1, sticky="W")

        # Manual input setup
        ttk.Label(frame, text="Manual Input:").grid(row=1, column=0, sticky="W", columnspan=2)

        # Polynomial selection
        ttk.Label(frame, text="Choose Polynomial:").grid(row=2, column=0, sticky="W")
        self.polynomial_var = tk.StringVar(value=list(self.get_polynomial_options().keys())[0])
        self.polynomial_menu = ttk.Combobox(frame, textvariable=self.polynomial_var, values=list(self.get_polynomial_options().keys()), state="readonly")
        self.polynomial_menu.grid(row=2, column=1, sticky="W")

        # Seed input
        ttk.Label(frame, text="Seed (e.g., 0b1010101):").grid(row=3, column=0, sticky="W")
        self.seed_entry = ttk.Entry(frame, width=50)
        self.seed_entry.grid(row=3, column=1, sticky="W")

        # Number of bits input
        ttk.Label(frame, text="Number of bits (n):").grid(row=4, column=0, sticky="W")
        self.n_entry = ttk.Entry(frame, width=50)
        self.n_entry.grid(row=4, column=1, sticky="W")

        # Buttons for running the simulation
        ttk.Button(frame, text="Run Simulation (Manual)", command=self.run_simulation).grid(row=5, column=0, pady=10, columnspan=2)

        # Output text area
        self.result_text = tk.Text(frame, width=80, height=10)
        self.result_text.grid(row=6, column=0, columnspan=2, pady=10)

        # Frame for graph
        self.graph_frame = ttk.Frame(self.root, padding="30")
        self.graph_frame.grid(row=1, column=0, sticky="NSEW")

    def get_polynomial_options(self):
        return {
            "f(x) = x^7 + x + 1": [1, 0, 0, 0, 0, 0, 1, 1],
            "f(x) = x^7 + x^5 + x^3 + 1": [1, 0, 1, 0, 1, 0, 0, 1],
            "f(x) = x^7 + x^6 + x^5 + x^2 + 1": [1, 1, 1, 0, 0, 1, 0, 1]
        }

    def run_simulation_from_file(self):
        try:
            # Read parameters from input.txt
            params = self.read_input_file("input.txt")
            polynomial_str = self.polynomial_var.get()
            seed_str = params["seed"]
            n = int(params["n"])

            # Get selected polynomial
            polynomial = self.get_polynomial_options().get(polynomial_str)
            if not polynomial:
                raise ValueError("Invalid polynomial selected.")

            # Validate seed
            if not seed_str.startswith("0b"):
                raise ValueError("Seed must be in binary format (e.g., 0b1010101).")
            seed = int(seed_str, 2)
            if seed <= 0 or seed >= (1 << len(polynomial)):
                raise ValueError(f"Seed must be a binary number with {len(polynomial)} bits.")

            # Run simulation
            sequence, period, S_star, p_value, hypothesis, critical_value, K = self.run_lfsr_simulation(polynomial, seed, n)

            # Write results to output.txt
            self.write_output_file("output.txt", polynomial_str, seed_str, n, sequence, period, S_star, p_value, hypothesis, critical_value, K)

            # Update text output
            self.update_result_text("output.txt", polynomial_str, seed_str, n, period, hypothesis)
            self.display_graph(sequence, period, S_star, p_value, hypothesis)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def run_simulation(self):
        try:
            polynomial_str = self.polynomial_var.get()
            seed_str = self.seed_entry.get()
            n = int(self.n_entry.get())

            polynomial = self.get_polynomial_options().get(polynomial_str)
            if not polynomial:
                raise ValueError("Invalid polynomial selected.")
            if not seed_str.startswith("0b"):
                raise ValueError("Seed must be in binary format (e.g., 0b1010101).")
            seed = int(seed_str, 2)
            if seed <= 0 or seed >= (1 << len(polynomial)):
                raise ValueError(f"Seed must be a binary number with {len(polynomial)} bits.")

            sequence, period, S_star, p_value, hypothesis, critical_value, K = self.run_lfsr_simulation(polynomial, seed, n)

            self.write_output_file("output.txt", polynomial_str, seed_str, n, sequence, period, S_star, p_value, hypothesis, critical_value, K)
            self.update_result_text("output.txt", polynomial_str, seed_str, n, period, hypothesis)
            self.display_graph(sequence, period, S_star, p_value, hypothesis)
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def read_input_file(self, filename):
        params = {}
        with open(filename, "r") as file:
            for line in file:
                key, value = line.strip().split("=")
                params[key.strip()] = value.strip()
        return params

    def write_output_file(self, filename, polynomial_str, seed_str, n, sequence, period, S_star, p_value, hypothesis, critical_value, K):
        with open(filename, "w") as f:
            f.write(f"Polynomial: {polynomial_str}\n")
            f.write(f"Seed: {seed_str}\n")
            f.write(f"Number of bits (n): {n}\n")
            f.write(f"Sequence: {''.join(map(str, sequence))}\n")
            f.write(f"Period: {period}\n")
            f.write(f"Chi-square statistic (S*): {S_star:.4f}\n")
            f.write(f"p-value: {p_value:.4f}\n")
            f.write(f"Critical value: {critical_value:.4f}\n")
            f.write(f"Number of intervals (K): {K}\n")
            f.write(f"Hypothesis test result: {hypothesis}\n")

    def update_result_text(self, filename, polynomial_str, seed_str, n, period, hypothesis):
        with open(filename, "r") as f:
            content = f.read()
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"Results saved to {filename}\n")
        self.result_text.insert(tk.END, content)

    def run_lfsr_simulation(self, polynomial, seed, n):
        state = seed
        sequence = []
        for _ in range(n):
            new_bit = sum(((state >> i) & 1) * polynomial[i] for i in range(len(polynomial))) % 2
            state = ((state >> 1) | (new_bit << (len(polynomial) - 1))) & ((1 << len(polynomial)) - 1)
            sequence.append(new_bit)
        period = self.calculate_period(sequence)
        S_star, p_value, hypothesis, critical_value, K = self.calculate_chi_square(sequence)
        return sequence, period, S_star, p_value, hypothesis, critical_value, K

    def calculate_period(self, sequence):
        n = len(sequence)
        for p in range(1, n // 2 + 1):  # Проверяем длину периода от 1 до n//2
            if n % p == 0:  # Только длины, которые могут быть делителями длины последовательности
                is_periodic = True
                for i in range(p, n):
                    if sequence[i] != sequence[i % p]:  # Проверяем, повторяется ли шаблон
                        is_periodic = False
                        break
                if is_periodic:
                    return p  # Возвращаем первый найденный период
        return n  # Если период не найден, возвращаем длину последовательности

    def calculate_chi_square(self, sequence):
        N = len(sequence)
        M = max(sequence) + 1
        K = M if M < 10 else int(5 * np.log10(N))
        intervals = np.linspace(0, M, K + 1)
        counts, _ = np.histogram(sequence, bins=intervals)
        v_i = counts
        P_i = 1 / K
        S_star = N * np.sum(((v_i / N - P_i) ** 2) / P_i)
        r = K - 1
        critical_value = chi2.ppf(0.95, df=r)
        p_value = 1 - chi2.cdf(S_star, df=r)
        hypothesis = "Rejected" if S_star > critical_value else "Not Rejected"
        return S_star, p_value, hypothesis, critical_value, K

    def display_graph(self, sequence, period, S_star, p_value, hypothesis):
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(sequence, drawstyle="steps-pre", marker="o", label="Sequence")
        ax.set_title(f"Sequence Visualization\nPeriod: {period} | S*: {S_star:.2f}, p-value: {p_value:.4f}\nHypothesis: {hypothesis}")
        ax.set_xlabel("Step")
        ax.set_ylabel("Bit Value")
        ax.legend()

        for widget in self.graph_frame.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.graph_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()


# Main execution
if __name__ == "__main__":
    root = tk.Tk()
    app = LFSRApp(root)
    root.mainloop()
