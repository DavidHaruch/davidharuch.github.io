import matplotlib.pyplot as plt
import numpy as np

# Set publication quality style for an academic thesis (Complete Sans-Serif)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica', 'Lucida Grande'],
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    # Force the math engine to use the exact same sans-serif fonts for numbers/symbols
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'DejaVu Sans',
    'mathtext.it': 'DejaVu Sans:italic',
    'mathtext.bf': 'DejaVu Sans:bold',
})

# ==========================================
# 1. INPUT PARAMETERS & SOLVER
# ==========================================
Q = 10            # Input Heat Dissipation/Wattage (W)
T_ambient = 25.0  # Ambient Temperature (noted as T3) (°C)

# Thermal Resistances (°C/W)
R1 = 3.7        
R2 = 0.5        
R3 = 1.0        

# Calculate backward from the ambient node (T3)
T3 = T_ambient
T2 = T3 + (Q * R3)
T1 = T2 + (Q * R2)
T0 = T1 + (Q * R1)

# ==========================================
# 2. SCHEMATIC GENERATION
# ==========================================
def draw_resistor(ax, x1, x2, y, h=0.12):
    """Draws a standard continuous schematic resistor from x1 to x2."""
    length = x2 - x1
    lead = 0.25 * length
    
    x_res = np.linspace(x1 + lead, x2 - lead, 7)
    y_res = [y, y + h, y - h, y + h, y - h, y + h, y]
    
    x_all = [x1, x1 + lead] + list(x_res[1:-1]) + [x2 - lead, x2]
    y_all = [y, y] + list(y_res[1:-1]) + [y, y]
    
    ax.plot(x_all, y_all, 'k-', lw=1.2)

# Initialize Figure
fig, ax = plt.subplots(figsize=(7, 2.2))

# Horizontal node coordinates along the chain
nodes = [1.0, 3.2, 5.4, 7.6]
y_val = 0.0

# Draw continuous resistors matching node intervals exactly
draw_resistor(ax, nodes[0], nodes[1], y_val)
draw_resistor(ax, nodes[1], nodes[2], y_val)
draw_resistor(ax, nodes[2], nodes[3], y_val)

# Plot standard open-terminal node circles
ax.plot(nodes, [y_val]*4, 'ko', markersize=5, markerfacecolor='white', markeredgewidth=1.2, zorder=5)

# Formal text placement directly below the respective nodes (No serif strings)
ax.text(nodes[0], y_val - 0.25, f'T$_0$ = {T0:.1f} °C', ha='center', va='top')
ax.text(nodes[1], y_val - 0.25, f'T$_1$ = {T1:.1f} °C', ha='center', va='top')
ax.text(nodes[2], y_val - 0.25, f'T$_2$ = {T2:.1f} °C', ha='center', va='top')
ax.text(nodes[3], y_val - 0.25, f'T$_3$ = {T3:.1f} °C', ha='center', va='top')

# Thermal resistance parameter labels placed right above each respective component
ax.text((nodes[0] + nodes[1])/2, y_val + 0.22, f'R$_1$ = {R1:.1f} °C/W', ha='center', va='bottom')
ax.text((nodes[1] + nodes[2])/2, y_val + 0.22, f'R$_2$ = {R2:.1f} °C/W', ha='center', va='bottom')
ax.text((nodes[2] + nodes[3])/2, y_val + 0.22, f'R$_3$ = {R3:.1f} °C/W', ha='center', va='bottom')

# Minimalist horizontal line with arrow indicating heat flux direction
ax.annotate('', xy=(nodes[3] + 0.4, y_val + 0.65), xytext=(nodes[0] - 0.4, y_val + 0.65),
            arrowprops=dict(arrowstyle="->", color="black", lw=1.0, shrinkA=0, shrinkB=0))
ax.text((nodes[0] + nodes[3])/2, y_val + 0.75, f'Q = {Q:.2f} W', ha='center', va='bottom')

# Clean layout boundaries
ax.set_xlim(nodes[0] - 1.0, nodes[3] + 1.0)
ax.set_ylim(-0.8, 1.3)
ax.axis('off')

plt.tight_layout()
plt.savefig('thermal_network_sans_serif.png', dpi=300, bbox_inches='tight')
# plt.close()
print("Figure successfully exported using uniform sans-serif configuration.")