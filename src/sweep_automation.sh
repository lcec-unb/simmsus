#!/usr/bin/env bash
set -euo pipefail

# Força locale "C" para garantir ponto decimal em printf/awk
export LC_ALL=C

echo "=============================================="
echo " SIMMSUS - Varredura 1 parâmetro (real)"
echo "=============================================="
echo ""

read -rp "Pasta BASE (contém simmsus.ex e simconfig.dat): " BASEDIR

# Normaliza caminho removendo barra final
BASEDIR="${BASEDIR%/}"

EXE="${BASEDIR}/simmsus.ex"
CFG="${BASEDIR}/simconfig.dat"

if [[ ! -d "$BASEDIR" ]]; then
  echo "ERRO: pasta não encontrada: $BASEDIR"
  exit 1
fi
if [[ ! -f "$EXE" ]]; then
  echo "ERRO: não encontrei $EXE"
  exit 1
fi
if [[ ! -f "$CFG" ]]; then
  echo "ERRO: não encontrei $CFG"
  exit 1
fi

echo ""
echo "Cole a CHAVE EXATA da linha (até os dois pontos), por exemplo:"
echo "  ALPHA .................................:"
echo "ou"
echo "  VOLUME FRACTION OF PARTICLES...........:"
echo ""
read -r -p "CHAVE: " KEY

if ! grep -Fq "$KEY" "$CFG"; then
  echo "ERRO: não encontrei essa chave no simconfig.dat: $KEY"
  exit 1
fi

echo ""
read -rp "Número de casos (inteiro >=2): " NCASE
read -rp "Valor mínimo (real): " VMIN
read -rp "Valor máximo (real): " VMAX

if [[ "$NCASE" -lt 2 ]]; then
  echo "ERRO: NCASE precisa ser >= 2"
  exit 1
fi

STAMP=$(date +"%Y%m%d_%H%M%S")
KEYNAME=$(echo "$KEY" | sed 's/[ .:]*$//' | awk '{print $1}')
ROOT="${BASEDIR}/sweep_${KEYNAME}_${STAMP}"

mkdir -p "$ROOT"

echo ""
echo ">> Pasta base: $BASEDIR"
echo ">> Criando varredura em: $ROOT"
echo ">> Parâmetro: $KEY"
echo ">> Casos: $NCASE  |  intervalo: [$VMIN, $VMAX]"
echo ""

fmt_sci () {
  # Ex: 1.0000E+01
  awk -v x="$1" 'BEGIN { printf("%.4E", x) }'
}

sanitize_name () {
  # Garante ponto decimal e remove caracteres problemáticos (espaços etc.)
  # Mantém apenas [A-Za-z0-9._+-Ee]
  echo "$1" | tr ',' '.' | tr -cd 'A-Za-z0-9._+-Ee'
}

for ((i=0; i<NCASE; i++)); do
  VAL=$(awk -v a="$VMIN" -v b="$VMAX" -v n="$NCASE" -v i="$i" \
        'BEGIN{ printf("%.16g", a + i*(b-a)/(n-1)) }')

  VALSCI=$(fmt_sci "$VAL")
  VALDIR=$(sanitize_name "$VALSCI")

  CASEDIR="${ROOT}/${KEYNAME}_${VALDIR}"
  mkdir -p "$CASEDIR"

  echo ">> Caso $((i+1))/$NCASE : $KEYNAME = $VALSCI"
  cp -f "$EXE" "$CASEDIR/simmsus.ex"
  cp -f "$CFG" "$CASEDIR/simconfig.dat"

  # Substitui apenas a parte após os dois pontos, preservando a chave
 awk -v key="$KEY" -v val="$VALSCI" '
    {
      lead = $0
      sub(/[^ \t].*$/, "", lead)   # prefixo de espaços
      line = $0
      sub(/^[ \t]+/, "", line)     # remove indentação para comparar
      if (index(line, key) == 1) {
        print lead key " " val
        next
      }
      print
    }
  ' "$CASEDIR/simconfig.dat" > "$CASEDIR/simconfig.dat.tmp" && mv "$CASEDIR/simconfig.dat.tmp" "$CASEDIR/simconfig.dat"


  # Execução desacoplada
  (
    cd "$CASEDIR"
    chmod +x simmsus.ex
    nohup ./simmsus.ex > run.log 2>&1 &
    echo "   PID: $!"
  )

done

echo ""
echo "=============================================="
echo "Varredura disparada!"
echo "Pasta raiz: $ROOT"
echo "Monitorar: tail -f $ROOT/<caso>/run.log"
echo "Listar PIDs: ps -fu \$USER | grep simmsus"
echo "=============================================="
