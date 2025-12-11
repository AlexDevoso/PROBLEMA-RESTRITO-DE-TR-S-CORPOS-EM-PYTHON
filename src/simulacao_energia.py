from funcoes_orbitais import *
from integradores import *



plt.rcParams['figure.dpi'] = 300
golden_ratio = (1 + np.sqrt(5)) / 2
figure_width = 3
figure_height = figure_width / golden_ratio
plt.rcParams['figure.figsize'] = [figure_width, figure_height]

plt.figure()
plt.plot(df_energia)
plt.xlabel("Passo de tempo")
plt.ylabel("Energia do sistema")
plt.title("Diferença de energia do sistema")
plt.grid(True)
plt.show()