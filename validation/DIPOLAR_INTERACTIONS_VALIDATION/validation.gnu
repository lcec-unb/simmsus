# ======================== validation.gnu ========================
# Executar: gnuplot validation.gnu
# Gera: comparacao.png
# Entrada esperada: magnetization_database.dat no formato Tecplot

infile  = "magnetization_database.dat"
outfile = "comparison.png"

# -------- Pré-processamento --------
# converte o arquivo Tecplot em blocos separados e extrai títulos das zonas
system("awk '/^[[:space:]]*[Vv]ariables/ {next} /^[[:space:]]*zone[[:space:]]+/ {if (printed>0) print \"\"; next} /^[[:space:]]*[0-9+\\.-]/ {print; printed=1; next} {next}' " . infile . " > .gp_blocks.dat")
system("awk '/^[[:space:]]*zone[[:space:]]+t=\"([^\"]+)\"/ {match($0,/t=\"([^\"]+)\"/,m); print m[1]}' " . infile . " > .gp_titles.txt")

# lê títulos dinamicamente
gettitle(i) = system(sprintf("sed -n '%dp' .gp_titles.txt | tr -d '\\n'", i))

T1 = gettitle(1); if (strlen(T1)==0) T1 = "SIMMSUS"
T2 = gettitle(2); if (strlen(T2)==0) T2 = "Langevin"
T3 = gettitle(3); if (strlen(T3)==0) T3 = "Ivanov - f^2"
T4 = gettitle(4); if (strlen(T4)==0) T4 = "Ivanov - f^3"

data = ".gp_blocks.dat"

# -------- Aparência --------
set term pngcairo size 1400,900 enhanced font "Helvetica,16"
set output outfile
set border lw 1.5
set grid back lw 1.0 dt 3
set key top left Left reverse spacing 1.2
set tics nomirror
set format x "%.2f"
set format y "%.4f"
set xlabel "Fraç\u00E3o volum\u00E9trica, {/Symbol f}" offset 0,0.5
set ylabel "M_{0}/M_{d}" offset 1.2,0

# estilos
set style line 1  lc rgb "#1f77b4" pt 7 ps 1.7 lw 2   # SIMMSUS (pontos+erro)
set style line 11 lc rgb "#1f77b4" lt 1 lw 2           # SIMMSUS (linha)
set style line 2  lc rgb "#d62728"  lt 1 lw 3          # Langevin
set style line 3  lc rgb "#2ca02c"  lt 1 lw 3          # Ivanov f^2
set style line 4  lc rgb "#ff7f0e"  lt 1 lw 3          # Ivanov f^3

# -------- Plot --------
plot \
    data index 0 using 1:2:3 with yerrorbars ls 1  title T1, \
    data index 0 using 1:2    with lines      ls 11 notitle, \
    data index 1 using 1:2    with lines      ls 2  title T2, \
    data index 2 using 1:2    with lines      ls 3  title T3, \
    data index 3 using 1:2    with lines      ls 4  title T4

unset output
print sprintf("Figura gerada: %s", outfile)
# ===============================================================
