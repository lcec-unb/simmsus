#!/usr/bin/env bash
# Analisa magnetization.plt em phi_0.01..phi_0.20
# Usa as últimas 20 amostras da 5ª coluna (Mz), multiplica por phi
# e grava "phi  media  desvio" em simmsus_magnetization.dat.

set -euo pipefail
export LC_ALL=C

OUT="simmsus_magnetization.dat"
: > "$OUT"

for p in $(seq 0.01 0.01 0.20); do
  phi=$(printf "%.2f" "$p")
  dir="phi_${phi}"
  file="${dir}/magnetization.plt"

  if [[ ! -f "$file" ]]; then
    echo "Aviso: arquivo não encontrado: $file" >&2
    # ainda assim escreve uma linha para manter 20 saídas
    printf "%.2f 0.0000E+00 0.0000E+00\n" "$phi" >> "$OUT"
    continue
  fi

  awk -v phi="$phi" '
    BEGIN {
      cap = 20; n = 0; idx = 0;
    }
    # pula cabeçalho do Tecplot
    /^ *[Vv]ariables *=/ { next }

    {
      # garante que temos ao menos 5 campos e que o 5º parece numérico
      # (aceita formatos 1.23E+02, 0.123, etc.)
      if (NF >= 5 && $5 ~ /^[ \t]*[+-]?[0-9]*\.?[0-9]+([Ee][+-]?[0-9]+)?[ \t]*$/) {
        # converte para número (awk lida bem com "E")
        val = $5 + 0
        # buffer circular das últimas 20 amostras
        buf[idx % cap] = val * phi
        idx++
        if (n < cap) n++
      }
    }
    END {
      if (n == 0) {
        printf("%.2f 0.0000E+00 0.0000E+00\n", phi)
        exit
      }
      # média e desvio-padrão amostral
      sum = 0.0; sumsq = 0.0
      for (i = 0; i < n; i++) {
        x = buf[i]
        sum   += x
        sumsq += x*x
      }
      mean = sum / n
      if (n > 1) { var = (sumsq - n*mean*mean) / (n - 1) } else { var = 0.0 }
      std = (var > 0 ? sqrt(var) : 0.0)

      printf("%.2f %.4E %.4E\n", phi, mean, std)
    }
  ' "$file" >> "$OUT"
done

echo "✅ Resultados salvos em $OUT"
