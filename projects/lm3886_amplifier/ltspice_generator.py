import os

def generate_ltspice_schematic(filename="current_loop.asc"):
    lines = []
    
    # 1. Header and Sheet Settings
    lines.append("VERSION 4")
    lines.append("SHEET 1 880 680")
    
    # 2. Component Symbols (SYMBOL name x y rotation)
    # R_ref and R_fb_pi (Inputs to PI controller)
    lines.append("SYMBOL res 160 144 R90")
    lines.append("WINDOW 0 0 56 VLeft 2")
    lines.append("WINDOW 3 0 32 VLeft 2")
    lines.append("SYMATTR Inst R_ref")
    lines.append("SYMATTR Value 10k")
    
    lines.append("SYMBOL res 160 240 R90")
    lines.append("WINDOW 0 0 56 VLeft 2")
    lines.append("WINDOW 3 0 32 VLeft 2")
    lines.append("SYMATTR Inst R_fb_pi")
    lines.append("SYMATTR Value 10k")

    # PI Op-Amp (OP27)
    lines.append("SYMBOL op27 320 192 R0")
    lines.append("SYMATTR Inst X_PI_OPAMP")
    
    # PI Feedback Components (R_z, R_p, C_i)
    lines.append("SYMBOL res 272 64 R90")
    lines.append("SYMATTR Inst R_z")
    lines.append("SYMATTR Value 10k")
    
    lines.append("SYMBOL res 208 16 R90")
    lines.append("SYMATTR Inst R_p")
    lines.append("SYMATTR Value 1k")
    
    lines.append("SYMBOL cap 336 0 R90")
    lines.append("SYMATTR Inst C_i")
    lines.append("SYMATTR Value 4.7n")

    # Voltage Divider (R_div1, R_div2)
    lines.append("SYMBOL res 496 176 R90")
    lines.append("SYMATTR Inst R_div1")
    lines.append("SYMATTR Value 4.7k")
    
    lines.append("SYMBOL res 544 240 R0")
    lines.append("SYMATTR Inst R_div2")
    lines.append("SYMATTR Value 1k")

    # Power Amp Placeholder (LM3886 using standard 5-pin opamp2 symbol)
    lines.append("SYMBOL opamp2 656 208 R0")
    lines.append("SYMATTR Inst X_POWER_AMP")
    lines.append("SYMATTR Value LM3886")
    
    # Power Amp Gain Network (R_pa1, R_pa2)
    lines.append("SYMBOL res 624 320 R0")
    lines.append("SYMATTR Inst R_pa1")
    lines.append("SYMATTR Value 1k")
    
    lines.append("SYMBOL res 688 352 R90")
    lines.append("SYMATTR Inst R_pa2")
    lines.append("SYMATTR Value 15k")

    # Load & Shunt (L_coil, R_coil, R_sh)
    lines.append("SYMBOL ind 800 224 R0")
    lines.append("SYMATTR Inst L_coil")
    lines.append("SYMATTR Value 0.1486m")
    
    lines.append("SYMBOL res 800 320 R0")
    lines.append("SYMATTR Inst R_coil_int")
    lines.append("SYMATTR Value 0.5")
    
    lines.append("SYMBOL res 800 416 R0")
    lines.append("SYMATTR Inst R_shunt")
    lines.append("SYMATTR Value 0.1")

    # Current Sense Op-Amp (OP27)
    lines.append("SYMBOL op27 656 528 R0")
    lines.append("SYMATTR Inst X_CS_OPAMP")
    
    # Current Sense Gain Network (R_cs1, R_cs2)
    lines.append("SYMBOL res 624 624 R0")
    lines.append("SYMATTR Inst R_cs1")
    lines.append("SYMATTR Value 1k")
    
    lines.append("SYMBOL res 688 656 R90")
    lines.append("SYMATTR Inst R_cs2")
    lines.append("SYMATTR Value 9k")

    # 3. Grounds (FLAG x y 0)
    lines.append("FLAG 272 240 0")     # PI non-inverting input GND
    lines.append("FLAG 544 320 0")     # Divider GND
    lines.append("FLAG 624 400 0")     # PA gain GND
    lines.append("FLAG 800 512 0")     # Shunt GND
    lines.append("FLAG 624 704 0")     # CS gain GND
    lines.append("FLAG 272 0 0")       # PI rails/GND helper if needed

    # 4. Wiring (WIRE x1 y1 x2 y2)
    # Reference input to R_ref
    lines.append("WIRE 48 160 96 160")
    lines.append("TEXT 32 160 Left 2!   ref")
    
    # R_ref and R_fb_pi summing into PI Inverting Input
    lines.append("WIRE 160 160 272 160")
    lines.append("WIRE 160 256 240 256")
    lines.append("WIRE 240 256 240 160")
    lines.append("WIRE 272 160 272 208")
    lines.append("WIRE 272 208 320 208")
    
    # PI Feedback loop routing
    lines.append("WIRE 240 160 240 80")
    lines.append("WIRE 240 80 208 80")
    lines.append("WIRE 208 32 208 16")
    lines.append("WIRE 208 16 288 16")
    lines.append("WIRE 336 16 336 0")
    lines.append("WIRE 336 0 432 0")
    lines.append("WIRE 432 0 432 224")
    lines.append("WIRE 400 224 432 224")
    
    # PI Output to Voltage Divider
    lines.append("WIRE 432 224 448 224")
    lines.append("WIRE 496 192 544 192")
    lines.append("WIRE 544 192 544 240")
    lines.append("WIRE 544 192 656 224") # To LM3886 +pin

    # LM3886 Output to Load Coil
    lines.append("WIRE 736 240 800 240")
    lines.append("WIRE 800 240 800 224")
    
    # Coil to Shunt and CS Amp Input
    lines.append("WIRE 800 304 800 320")
    lines.append("WIRE 800 400 800 416")
    lines.append("WIRE 800 400 592 400")
    lines.append("WIRE 592 400 592 544")
    lines.append("WIRE 592 544 656 544") # To CS Amp +pin

    # Long Feedback path from CS Amp Output back to PI Input
    lines.append("WIRE 736 560 768 560")
    lines.append("WIRE 768 560 768 480")
    lines.append("WIRE 768 480 96 480")
    lines.append("WIRE 96 480 96 256")
    lines.append("WIRE 96 256 112 256")

    # 5. Simulation Directives
    lines.append("TEXT -32 -96 Left 2! .trans 1n 500u 0 10n")
    
    # Write to file
    with open(filename, "w") as f:
        f.write("\n".join(lines))
    print(f"Successfully generated {filename}")

if __name__ == "__main__":
    generate_ltspice_schematic()