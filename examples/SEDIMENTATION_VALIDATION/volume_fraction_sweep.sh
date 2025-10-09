#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║ 💻  AESC – Varredura de fração volumétrica (SIMMSUS)                         ║
# ║ 🛠️  Cria pastas, copia simconfig.dat/simmsus.ex, edita VOLUME FRACTION...     ║
# ║ 🚀  Dispara as simulações com nohup                                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

set -euo pipefail
export LC_ALL=C
export LC_NUMERIC=C

# -------- Configurações rápidas --------
# Estilo de nome de pasta (seguro para shell). Se quiser "φ=0.01", troque a linha abaixo.
NAME_STYLE="phi_%.2f"
BASE_DIR="$(pwd)"
SRC_CONFIG="simconfig.dat"
SRC_BIN="simmsus.ex"   # <- corrigido: dois 'm'

# -------- Verificações de pré-execução --------
if [[ ! -f "$SRC_CONFIG" ]]; then
  echo "ERRO: Arquivo '$SRC_CONFIG' não encontrado em: $BASE_DIR"
  exit 1
fi
if [[ ! -f "$SRC_BIN" ]]; then
  echo "ERRO: Arquivo '$SRC_BIN' não encontrado em: $BASE_DIR"
  exit 1
fi
chmod +x "$SRC_BIN" || true

echo ""
echo ">>> Iniciando varredura de φ de 0.01 a 0.20 (passo 0.01) ..."
echo ""

# -------- Loop principal --------
for n in $(seq 1 20); do
  # φ como decimal e em notação científica com 4 casas decimais (ex.: 1.0000E-02)
  phi_dec=$(awk -v k="$n" 'BEGIN{printf("%.2f", k/100.0)}')
  phi_fmt=$(awk -v x="$phi_dec" 'BEGIN{printf("%.4E", x)}')

  # Nome da pasta
  folder=$(printf "$NAME_STYLE" "$phi_dec")

  # Criação da pasta
  mkdir -p "$folder"

  # Copiar arquivos
  cp -f "$SRC_CONFIG" "$folder/"
  cp -f "$SRC_BIN"    "$folder/"

  # Editar linha "VOLUME FRACTION OF PARTICLES...........: <valor>"
  # Mantém exatamente o prefixo (com pontos e dois pontos) e substitui apenas o valor.
  # Ex.: " VOLUME FRACTION OF PARTICLES...........: 1.5000E-01"
  sed -E -i \
    's/^([[:space:]]*VOLUME[[:space:]]+FRACTION[[:space:]]+OF[[:space:]]+PARTICLES[.]*:[[:space:]]*)([+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?)/\1'"$phi_fmt"'/' \
    "$folder/$SRC_CONFIG"

  # Disparar simulação com nohup; logs e PID por pasta
  (
    cd "$folder"
    nohup ./simmsus.ex > "log.simmsus" 2>&1 &
    pid=$!
    echo "$pid" > "simmsus.pid"
    echo "✓ Pasta '$folder' | VOLUME FRACTION OF PARTICLES = $phi_fmt | PID $pid"
  )
done

echo ""
echo ">>> Varredura iniciada. Todas as simulações foram disparadas com nohup."
echo "    Acompanhe com: tail -f */log.simmsus   ou veja PIDs em: */simmsus.pid"
