#!/bin/bash
# ======================================================================
# SIMMSUS - Terminal configurator (pure bash UI) + simconfig.dat writer
# ======================================================================

set -Eeuo pipefail

# ========= único helper de UI (como combinado) ========================
print_header() {
  clear
  echo ""
  echo "╔══════════════════════════════════════════════════════════════════════════════╗"
  echo "║             🧲 SIMMSUS - SIMULATION OF MAGNETIC SUSPENSIONS 🧲               ║"
  echo "╠══════════════════════════════════════════════════════════════════════════════╣"
  echo "║                      💻 simconfig.dat generator 💻                           ║"
  echo "╚══════════════════════════════════════════════════════════════════════════════╝"
  echo ""
}

# ========= helpers de formatação (para o writer) ==============================
# (não são UI; apenas garantem o formato do arquivo)
sci () {
  local val="${1:-0}"
  val="${val//,/\.}"                  # vírgula -> ponto
  LC_NUMERIC=C printf '%.4E' "$val" | tr 'e' 'E'
}
zpad () {
  local w="$1"; shift
  local v="${1:-0}"
  # garante que é inteiro (remove qualquer coisa não-numérica, exceto '-')
  v="${v//[^0-9-]/}"
  [[ "$v" =~ ^-?[0-9]+$ ]] || v=0
  LC_NUMERIC=C printf "%0${w}d" "$v"
}
oneline () { printf '%s' "${1:-}" | tr -d '\r\n'; }

bool_from_choice () {  # 0 -> TRUE, 1 -> FALSE
  case "${1:-1}" in
    0) echo TRUE  ;;
    1) echo FALSE ;;
    *) echo FALSE ;;
  esac
}

# ========= inicialização de variáveis (evita "unbound") ======================
choice1=""; lambda=""
choice1_0=""; choice1_05=""; choice1_06=""; continuefrom=""
choice1_07=""
number_particles=""; number_rea=""; sim_time=""; time_step=""; write_interval=""
phi=""; aspect_ratio=""
choice_rec_pos=""; choice_rec_vel=""; choice_rec_dip=""; choice_ovito=""
percentnonmagpart=""
choice_ordered_arr=""; choice_init_sph=""; choice_ordered_dip=""
choice1_1=""; choice1_2=""; choice1_3=""
choice2=""; pecletnumber=""
choice3=""; stokesnumber=""; choice3_2=""
choice4=""
choice5=""; choice6=""; shearamplitude=""; shearfrequency=""
choice7=""; choice8=""
alpha=""; omega_field=""; omega_field2=""; c1_duffing=""; c2_duffing=""; c3_duffing=""; c4_duffing=""
max_omega=""; numberintdyn=""

# ==========================
# ========== MENUS =========
# ==========================

# ===== MENU 1: particle type ===================================================
# (cores simples; sem funções novas)
CLR_NC="\033[0m"     # reset
CLR_TITLE="\033[1;36m"
CLR_LABEL="\033[1;37m"
CLR_DIM="\033[2m"
CLR_ERR="\033[1;31m"

while :; do
  print_header
   echo -e "${CLR_TITLE}Particle data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Choose one of the following options:${CLR_NC}"
  echo
  echo "  ┌───────────────────────────────┐"
  echo "  │  1) Magnetic particles        │"
  echo "  │  2) Non-magnetic particles    │"
  echo "  │  3) A mixture of both         │"
  echo "  └───────────────────────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}Select [1-3]:${CLR_NC} ) " choice1
  case "$choice1" in
    1|2|3) break ;;
    *) echo -e "${CLR_ERR}Invalid choice. Please type 1, 2 or 3.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 1.01: lambda if magnetic or mixture ===============================
case "$choice1" in
  1|3)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Magnetic effects${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Dipolar coupling parameter${CLR_NC}"
      echo
      echo "  ┌───────────────────────────────┐"
      echo "  │   Enter the value for λ       │"
      echo "  │         (λ or Ψ_m)            │"
      echo "  └───────────────────────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}λ or Ψ_m:${CLR_NC} ) " lambda
      # Validação simples (numérico)
      if [[ $lambda =~ ^[0-9]*\.?[0-9]+$ ]]; then
        break
      else
        echo -e "${CLR_ERR}Invalid value. Please enter a number.${CLR_NC}"
        sleep 1
      fi
    done
    clear
  ;;
esac

# ===== MENU 1.0: fluctuation mode =============================================
while :; do
  print_header
  echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Non-inertial observer (fluctuation mode)?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice1_0=0; break ;;
    N|NO)  choice1_0=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 1.05: Monte Carlo (drives STATIC vs DYNAMIC) ======================
while :; do
  print_header
  echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Monte Carlo simulations?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice1_05=0; break ;;
    N|NO)  choice1_05=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 1.06: continue old simulation? ====================================
while :; do
  print_header
  echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Do you want to continue an old simulation?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice1_06=0; break ;;
    N|NO)  choice1_06=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 1.06.a: "Continue from iteration" (only if YES above) =============
if [ "${choice1_06:-1}" -eq 0 ]; then
  while :; do
    print_header
    echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
    echo
    echo -e "${CLR_LABEL}Continue simulation data${CLR_NC}"
    echo
    echo "  ┌───────────────────────────────────────────────┐"
    echo "  │  Enter the iteration number to resume from    │"
    echo "  │  (must be an integer ≥ 0)                     │"
    echo "  └───────────────────────────────────────────────┘"
    echo
    read -r -p "$(echo -e ${CLR_DIM}Continue from iteration:${CLR_NC} ) " continuefrom
    if [[ "$continuefrom" =~ ^[0-9]+$ ]]; then
      break
    else
      echo -e "${CLR_ERR}Invalid value. Please enter a non-negative integer.${CLR_NC}"
      sleep 1
    fi
  done
  clear
fi

# ===== MENU 1.07: polydispersity? =============================================
while :; do
  print_header
  echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Do you want to consider particle polydispersity?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice1_07=0; break ;;
    N|NO)  choice1_07=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 1.1: essential data (multi-field form) ============================
while :; do
  print_header
  echo -e "${CLR_TITLE}Essential data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Simulation parameters${CLR_NC}"
  echo
  echo "  ┌───────────────────────────────────────────────┐"
  echo "  │          Fill in the fields below             │"
  echo "  └───────────────────────────────────────────────┘"
  echo

  # Number of particles (int)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Number of particles:${CLR_NC} ) " number_particles
    if [[ "$number_particles" =~ ^[0-9]+$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a non-negative integer.${CLR_NC}"; sleep 1
  done

  # Number of experiments (int)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Number of experiments:${CLR_NC} ) " number_rea
    if [[ "$number_rea" =~ ^[0-9]+$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a non-negative integer.${CLR_NC}"; sleep 1
  done

  # Simulation time (float)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Simulation time:${CLR_NC} ) " sim_time
    sim_time="${sim_time//,/.}"
    if [[ "$sim_time" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a number (e.g., 10, 3.5, 1e-2).${CLR_NC}"; sleep 1
  done

  # Numerical time step (float)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Numerical time step:${CLR_NC} ) " time_step
    time_step="${time_step//,/.}"
    if [[ "$time_step" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a number (e.g., 0.01, 1e-3).${CLR_NC}"; sleep 1
  done

  # Write interval (int)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Write interval:${CLR_NC} ) " write_interval
    if [[ "$write_interval" =~ ^[0-9]+$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a non-negative integer.${CLR_NC}"; sleep 1
  done

  # Volume fraction (float)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Volume fraction:${CLR_NC} ) " phi
    phi="${phi//,/.}"
    if [[ "$phi" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a number (e.g., 0.05, 5e-2).${CLR_NC}"; sleep 1
  done

  # Box aspect ratio (int)
  while :; do
    read -r -p "$(echo -e ${CLR_DIM}Box aspect ratio:${CLR_NC} ) " aspect_ratio
    if [[ "$aspect_ratio" =~ ^[0-9]+$ ]]; then break; fi
    echo -e "${CLR_ERR}Enter a non-negative integer.${CLR_NC}"; sleep 1
  done

  break
done
clear

# ===== Output/recording options ===============================================
# Record positions
while :; do
  print_header
  echo -e "${CLR_TITLE}Output options${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Record particle positions to file?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_rec_pos=0; break ;;
    N|NO)  choice_rec_pos=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# Record velocities
while :; do
  print_header
  echo -e "${CLR_TITLE}Output options${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Record particle velocities to file?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_rec_vel=0; break ;;
    N|NO)  choice_rec_vel=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# Record dipoles
while :; do
  print_header
  echo -e "${CLR_TITLE}Output options${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Record particle dipoles to file?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_rec_dip=0; break ;;
    N|NO)  choice_rec_dip=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# OVITO vs Tecplot
while :; do
  print_header
  echo -e "${CLR_TITLE}Output options${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Use OVITO output format? (No = Tecplot)${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_ovito=0; break ;;
    N|NO)  choice_ovito=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 1.15: percent non-magnetic (for mixture) ==========================
case "$choice1" in
  3)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Particle distribution${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Percentage of non-magnetic particles${CLR_NC}"
      echo
      echo "  ┌───────────────────────────────────────────────┐"
      echo "  │  Enter fraction as a number (e.g., 0.50)      │"
      echo "  │  Meaning: 0.50 = 50% non-magnetic             │"
      echo "  └───────────────────────────────────────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}% of non-magnetic particles:${CLR_NC} ) " percentnonmagpart
      percentnonmagpart="${percentnonmagpart//,/.}"
      if [[ "$percentnonmagpart" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      else
        echo -e "${CLR_ERR}Invalid value. Please enter a number (e.g., 0.5).${CLR_NC}"
        sleep 1
      fi
    done
    clear
  ;;
esac

# ===== Initial configuration flags ============================================
# Ordered arrangement?
while :; do
  print_header
  echo -e "${CLR_TITLE}Initial configuration${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Ordered (T) or random (F) arrangement?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Ordered   │"
  echo "  │  n) Random    │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_ordered_arr=0; break ;;
    N|NO)  choice_ordered_arr=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# Initial spherical aggregate?
while :; do
  print_header
  echo -e "${CLR_TITLE}Initial configuration${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Initial spherical aggregate?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_init_sph=0; break ;;
    N|NO)  choice_init_sph=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# Ordered dipoles?
while :; do
  print_header
  echo -e "${CLR_TITLE}Initial configuration${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Ordered dipoles?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice_ordered_dip=0; break ;;
    N|NO)  choice_ordered_dip=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== Rotating particles? (magnetic/mixture) =================================
case "$choice1" in
  1|3)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Do you want to consider rotating particles?${CLR_NC}"
      echo
      echo "  ┌───────────────┐"
      echo "  │  y) Yes       │"
      echo "  │  n) No        │"
      echo "  └───────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
      case "${ans^^}" in
        Y|YES) choice1_1=0; break ;;
        N|NO)  choice1_1=1; break ;;
        *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== Periodic forces? (if magnetic/mixture) =================================
case "$choice1" in
  1|3)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Magnetic interactions${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Compute dipolar forces periodically?${CLR_NC}"
      echo
      echo "  ┌───────────────┐"
      echo "  │  y) Yes       │"
      echo "  │  n) No        │"
      echo "  └───────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
      case "${ans^^}" in
        Y|YES) choice1_2=0; break ;;
        N|NO)  choice1_2=1; break ;;
        *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== Periodic torques? (only if rotating==YES) ==============================
case "${choice1_1:-1}" in
  0)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Magnetic interactions${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Compute dipolar torques periodically?${CLR_NC}"
      echo
      echo "  ┌───────────────┐"
      echo "  │  y) Yes       │"
      echo "  │  n) No        │"
      echo "  └───────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
      case "${ans^^}" in
        Y|YES) choice1_3=0; break ;;
        N|NO)  choice1_3=1; break ;;
        *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== MENU 2: brownian? ======================================================
while :; do
  print_header
  echo -e "${CLR_TITLE}Particle data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Choose one of the following options:${CLR_NC}"
  echo
  echo "  1) Brownian particles"
  echo "  2) Non-Brownian particles"
  echo
  read -r -p "$(echo -e ${CLR_DIM}Select [1-2]:${CLR_NC} ) " choice2
  case "$choice2" in
    1|2) break ;;
    *) echo -e "${CLR_ERR}Please enter 1 or 2.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

case "$choice2" in
  1)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Brownian effects${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Péclet number of the particles${CLR_NC}"
      echo
      read -r -p "$(echo -e ${CLR_DIM}Pe:${CLR_NC} ) " pecletnumber
      pecletnumber="${pecletnumber//,/.}"
      if [[ "$pecletnumber" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      else
        echo -e "${CLR_ERR}Invalid value. Please enter a number (e.g., 12.5).${CLR_NC}"
        sleep 1
      fi
    done
    clear
  ;;
esac

# ===== MENU 3: inertia (massive vs non-massive) ===============================
while :; do
  print_header
  echo -e "${CLR_TITLE}Particle data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Choose one of the following options:${CLR_NC}"
  echo
  echo "  1) Massive particles"
  echo "  2) Non-massive particles"
  echo
  read -r -p "$(echo -e ${CLR_DIM}Select [1-2]:${CLR_NC} ) " choice3
  case "$choice3" in
    1|2) break ;;
    *) echo -e "${CLR_ERR}Please enter 1 or 2.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

case "$choice3" in
  1)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Inertial effects${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Particle Stokes number${CLR_NC}"
      echo
      read -r -p "$(echo -e ${CLR_DIM}St:${CLR_NC} ) " stokesnumber
      stokesnumber="${stokesnumber//,/.}"
      if [[ "$stokesnumber" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      else
        echo -e "${CLR_ERR}Invalid value. Please enter a number (e.g., 0.25).${CLR_NC}"
        sleep 1
      fi
    done
    clear
  ;;
  2)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Hydrodynamic interactions${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Do you want to consider hydrodynamic interactions?${CLR_NC}"
      echo
      echo "  ┌───────────────┐"
      echo "  │  y) Yes       │"
      echo "  │  n) No        │"
      echo "  └───────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
      case "${ans^^}" in
        Y|YES) choice3_2=0; break ;;
        N|NO)  choice3_2=1; break ;;
        *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== MENU 4: sedimentation? =================================================
while :; do
  print_header
  echo -e "${CLR_TITLE}Particle data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Choose one of the following options:${CLR_NC}"
  echo
  echo "  1) Neutrally buoyant particles"
  echo "  2) Sedimenting particles"
  echo
  read -r -p "$(echo -e ${CLR_DIM}Select [1-2]:${CLR_NC} ) " choice4
  case "$choice4" in
    1|2) break ;;
    *) echo -e "${CLR_ERR}Please enter 1 or 2.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 5: shear? =========================================================
while :; do
  print_header
  echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
  echo
  echo -e "${CLR_LABEL}Suspension subjected to shear-rate?${CLR_NC}"
  echo
  echo "  ┌───────────────┐"
  echo "  │  y) Yes       │"
  echo "  │  n) No        │"
  echo "  └───────────────┘"
  echo
  read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
  case "${ans^^}" in
    Y|YES) choice5=0; break ;;
    N|NO)  choice5=1; break ;;
    *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
  esac
done
clear

# ===== MENU 6: oscillatory shear? (if shear YES) ==============================
case "$choice5" in
  0)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Oscillatory shear imposed?${CLR_NC}"
      echo
      echo "  ┌───────────────┐"
      echo "  │  y) Yes       │"
      echo "  │  n) No        │"
      echo "  └───────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
      case "${ans^^}" in
        Y|YES) choice6=0; break ;;
        N|NO)  choice6=1; break ;;
        *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== MENUS 7/8: shear data ==================================================
if [ "${choice5:-1}" = "0" ]; then
case "${choice6:-1}" in
  0)
    # Amplitude + frequency
    while :; do
      print_header
      echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Shear rate data${CLR_NC}"
      echo
      read -r -p "$(echo -e ${CLR_DIM}Shear rate amplitude:${CLR_NC} ) " shearamplitude
      shearamplitude="${shearamplitude//,/.}"
      if [[ ! "$shearamplitude" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        echo -e "${CLR_ERR}Invalid value. Enter a number (e.g., 0.5).${CLR_NC}"
        sleep 1
        continue
      fi

      read -r -p "$(echo -e ${CLR_DIM}Shear rate frequency:${CLR_NC} ) " shearfrequency
      shearfrequency="${shearfrequency//,/.}"
      if [[ ! "$shearfrequency" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        echo -e "${CLR_ERR}Invalid value. Enter a number (e.g., 1.0).${CLR_NC}"
        sleep 1
        continue
      fi
      break
    done
    clear
  ;;
  1)
    # Only amplitude
    while :; do
      print_header
      echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Shear rate data${CLR_NC}"
      echo
      read -r -p "$(echo -e ${CLR_DIM}Shear rate amplitude:${CLR_NC} ) " shearamplitude
      shearamplitude="${shearamplitude//,/.}"
      if [[ "$shearamplitude" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      else
        echo -e "${CLR_ERR}Invalid value. Enter a number (e.g., 0.5).${CLR_NC}"
        sleep 1
      fi
    done
    clear
  ;;
esac
fi

# ===== MENU 9: external magnetic field (for magnetic/mixture) =================
case "$choice1" in
  1|3)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Suspension data${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}External magnetic field applied?${CLR_NC}"
      echo
      echo "  ┌───────────────┐"
      echo "  │  y) Yes       │"
      echo "  │  n) No        │"
      echo "  └───────────────┘"
      echo
      read -r -p "$(echo -e ${CLR_DIM}[y/n]:${CLR_NC} ) " ans
      case "${ans^^}" in
        Y|YES) choice7=0; break ;;
        N|NO)  choice7=1; break ;;
        *) echo -e "${CLR_ERR}Please type y or n.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== MENU 10: field type ====================================================
case "${choice7:-1}" in
  0)
    while :; do
      print_header
      echo -e "${CLR_TITLE}External field data${CLR_NC}"
      echo
      echo -e "${CLR_LABEL}Choose one of the following options for the applied field:${CLR_NC}"
      echo
      echo "  1) Steady-state field"
      echo "  2) Oscillatory field"
      echo "  3) Rotative field"
      echo "  4) Oscillatory field with double frequency"
      echo "  5) Nonlinear Duffing solution applied field"
      echo "  6) Field-ramp (Dynamical increase)"
      echo
      read -r -p "$(echo -e ${CLR_DIM}Select [1-6]:${CLR_NC} ) " choice8
      case "$choice8" in
        1|2|3|4|5|6) break ;;
        *) echo -e "${CLR_ERR}Please enter a number between 1 and 6.${CLR_NC}"; sleep 1 ;;
      esac
    done
    clear
  ;;
esac

# ===== MENUS 11.15: field parameters ========================================
choice8="${choice8//$'\r'/}" 
choice8="${choice8//[[:space:]]/}" 

choice8="$(printf '%s' "${choice8:-}" | tr -d '\r\n[:space:]')"

case "$choice8" in
  1)  # Steady-state
    while :; do
      print_header
      echo -e "${CLR_TITLE}Field data${CLR_NC}"
      echo
      printf "%b" "${CLR_DIM}Nondimensional field (alpha):${CLR_NC} "
      read -r alpha
      alpha="${alpha//,/.}"
      [[ "$alpha" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && break
      echo -e "${CLR_ERR}Invalid number. Example: 2.0${CLR_NC}"; sleep 1
    done
    clear
  ;;

  2)  # Oscillatory
    while :; do
      print_header
      echo -e "${CLR_TITLE}Field data${CLR_NC}"
      echo
      printf "%b" "${CLR_DIM}Nondimensional field (alpha):${CLR_NC} "
      read -r alpha; alpha="${alpha//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional field frequency:${CLR_NC} "
      read -r omega_field; omega_field="${omega_field//,/.}"
      if [[ "$alpha" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$omega_field" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      fi
      echo -e "${CLR_ERR}Invalid values. Examples: alpha=2.0, freq=0.2${CLR_NC}"; sleep 1
    done
    clear
  ;;

  3)  # Rotative
    while :; do
      print_header
      echo -e "${CLR_TITLE}Field data${CLR_NC}"
      echo
      printf "%b" "${CLR_DIM}Nondimensional field (alpha):${CLR_NC} "
      read -r alpha; alpha="${alpha//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional field frequency:${CLR_NC} "
      read -r omega_field; omega_field="${omega_field//,/.}"
      if [[ "$alpha" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$omega_field" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      fi
      echo -e "${CLR_ERR}Invalid values. Examples: alpha=2.0, freq=0.2${CLR_NC}"; sleep 1
    done
    clear
  ;;

  4)  # Double frequency
    while :; do
      print_header
      echo -e "${CLR_TITLE}Field data${CLR_NC}"
      echo
      printf "%b" "${CLR_DIM}Nondimensional field (alpha):${CLR_NC} "
      read -r alpha; alpha="${alpha//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional field frequency 1:${CLR_NC} "
      read -r omega_field; omega_field="${omega_field//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional field frequency 2:${CLR_NC} "
      read -r omega_field2; omega_field2="${omega_field2//,/.}"
      if [[ "$alpha" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$omega_field" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$omega_field2" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      fi
      echo -e "${CLR_ERR}Invalid values. Example: alpha=2.0, f1=0.2, f2=2.0${CLR_NC}"; sleep 1
    done
    clear
  ;;

  5)  # Duffing
    while :; do
      print_header
      echo -e "${CLR_TITLE}Field data — Duffing${CLR_NC}"
      echo
      printf "%b" "${CLR_DIM}Nondimensional field (alpha):${CLR_NC} "
      read -r alpha; alpha="${alpha//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional field frequency:${CLR_NC} "
      read -r omega_field; omega_field="${omega_field//,/.}"
      printf "%b" "${CLR_DIM}C1 Duffing coefficient:${CLR_NC} "
      read -r c1_duffing; c1_duffing="${c1_duffing//,/.}"
      printf "%b" "${CLR_DIM}C2 Duffing coefficient:${CLR_NC} "
      read -r c2_duffing; c2_duffing="${c2_duffing//,/.}"
      printf "%b" "${CLR_DIM}C3 Duffing coefficient:${CLR_NC} "
      read -r c3_duffing; c3_duffing="${c3_duffing//,/.}"
      printf "%b" "${CLR_DIM}C4 Duffing coefficient:${CLR_NC} "
      read -r c4_duffing; c4_duffing="${c4_duffing//,/.}"
      if [[ "$alpha" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$omega_field" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$c1_duffing" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$c2_duffing" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$c3_duffing" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$c4_duffing" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]; then
        break
      fi
      echo -e "${CLR_ERR}Invalid values. Use numeric inputs (e.g., 0.1).${CLR_NC}"; sleep 1
    done
    clear
  ;;

  6)  # Dynamic increase (ramp)
    while :; do
      print_header
      echo -e "${CLR_TITLE}Field data — Frequency ramp${CLR_NC}"
      echo
      printf "%b" "${CLR_DIM}Nondimensional field (alpha):${CLR_NC} "
      read -r alpha; alpha="${alpha//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional minimum frequency:${CLR_NC} "
      read -r omega_field; omega_field="${omega_field//,/.}"
      printf "%b" "${CLR_DIM}Nondimensional maximum frequency:${CLR_NC} "
      read -r max_omega; max_omega="${max_omega//,/.}"
      printf "%b" "${CLR_DIM}Number of intervals:${CLR_NC} "
      read -r numberintdyn
      if [[ "$alpha" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$omega_field" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$max_omega" =~ ^-?([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] && \
         [[ "$numberintdyn" =~ ^[0-9]+$ ]]; then
        break
      fi
      echo -e "${CLR_ERR}Invalid values. alpha/frequencies must be numeric; intervals must be integer.${CLR_NC}"
      sleep 1
    done
    clear
  ;;

  *) : ;;  # nada selecionado / sem campo externo
esac

# ============================================================
# SANITIZAÇÃO FINAL – GARANTE VALORES UNITÁRIOS EM VARIÁVEIS 
# NÃO UTILIZADAS NO FLUXO PARA EVITAR NAN NO SIMMSUS
# ============================================================

# Todas as variáveis numéricas que podem ficar vazias recebem valor 1
sanitize_num() {
    local v="$1"
    # se vazio, nulo ou apenas espaço → retorna 1
    if [[ -z "${v// }" ]]; then
        echo "1"
    else
        echo "$v"
    fi
}

# Todas as variáveis relevantes passam por sanitize_num
lambda=$(sanitize_num "$lambda")
pecletnumber=$(sanitize_num "$pecletnumber")
stokesnumber=$(sanitize_num "$stokesnumber")
shearamplitude=$(sanitize_num "$shearamplitude")
shearfrequency=$(sanitize_num "$shearfrequency")

alpha=$(sanitize_num "$alpha")
omega_field=$(sanitize_num "$omega_field")
omega_field2=$(sanitize_num "$omega_field2")
c1_duffing=$(sanitize_num "$c1_duffing")
c2_duffing=$(sanitize_num "$c2_duffing")
c3_duffing=$(sanitize_num "$c3_duffing")
c4_duffing=$(sanitize_num "$c4_duffing")
max_omega=$(sanitize_num "$max_omega")
numberintdyn=$(sanitize_num "$numberintdyn")

percentnonmagpart=$(sanitize_num "$percentnonmagpart")

# Quantidades e parâmetros numéricos
number_particles=$(sanitize_num "$number_particles")
number_rea=$(sanitize_num "$number_rea")
sim_time=$(sanitize_num "$sim_time")
time_step=$(sanitize_num "$time_step")
write_interval=$(sanitize_num "$write_interval")
phi=$(sanitize_num "$phi")
aspect_ratio=$(sanitize_num "$aspect_ratio")

# Parâmetros físicos derivados
continuefrom=$(sanitize_num "$continuefrom")


# ==============================
# ========== WRITER ============
# ==============================

# --- Deriva flags lógicas a partir das respostas ------------------------------
case "$choice1" in
  1|3) MAG_PART=TRUE ;; *) MAG_PART=FALSE ;; esac
case "$choice1" in
  3) MIX_TYPES=TRUE ;; *) MIX_TYPES=FALSE ;; esac

SOLVE_TORQUE=$(bool_from_choice "${choice1_1:-1}")
PERIODIC_FORCES=$(bool_from_choice "${choice1_2:-1}")
PERIODIC_TORQUES=$( [ "${choice1_1:-1}" = "0" ] && bool_from_choice "${choice1_3:-1}" || echo FALSE )

FLUCT_MODE=$(bool_from_choice "${choice1_0:-1}")

CONTINUE_OLD=$(bool_from_choice "${choice1_06:-1}")
CONT_FROM="${continuefrom:-10000}"

POLY=$(bool_from_choice "${choice1_07:-1}")

if [ "${choice2:-2}" = "1" ]; then BROWNIAN=TRUE; else BROWNIAN=FALSE; fi

if [ "${choice3:-2}" = "1" ]; then PARTICLE_INERTIA=TRUE; else PARTICLE_INERTIA=FALSE; fi

if [ "${choice3:-2}" = "2" ]; then
  HDI=$(bool_from_choice "${choice3_2:-0}")
else
  HDI=TRUE
fi

if [ "${choice4:-1}" = "2" ]; then SEDIMENTATION=TRUE; else SEDIMENTATION=FALSE; fi

TURN_ON_SHEAR=$(bool_from_choice "${choice5:-1}")
OSC_SHEAR=$( [ "${choice5:-1}" = "0" ] && bool_from_choice "${choice6:-1}" || echo FALSE )
DYN_SHEAR=FALSE
GAMMA_D="${shearamplitude:-0}"
FREQ_SHEAR="${shearfrequency:-0}"

APPLY_FIELD=$(bool_from_choice "${choice7:-1}")
FIELD_POS=2

OSC_FIELD=FALSE; ROT_FIELD=FALSE; DOUBLE_FREQ=FALSE; DUFFING=FALSE; DYN_FREQ=FALSE
alpha_field="${alpha:-0}"
omega1="${omega_field:-0}"
omega2="${omega_field2:-0}"
c1="${c1_duffing:-0}"; c2="${c2_duffing:-0}"; c3="${c3_duffing:-0}"; c4="${c4_duffing:-0}"
omega_max="${max_omega:-0}"; n_intervals="${numberintdyn:-0}"

if [ "$APPLY_FIELD" = "TRUE" ]; then
  case "${choice8:-}" in
    2) OSC_FIELD=TRUE ;;
    3) ROT_FIELD=TRUE ;;
    4) DOUBLE_FREQ=TRUE ;;
    5) DUFFING=TRUE ;;
    6) DYN_FREQ=TRUE ;;
  esac
else
  alpha_field=0; omega1=0; omega2=0; c1=0; c2=0; c3=0; c4=0; omega_max=0; n_intervals=0
fi

# Quantidades e parâmetros numéricos
N_PART="${number_particles:-500}"
N_REALZ="${number_rea:-1}"
T_SIM="${sim_time:-1e+1}"
DT="${time_step:-1e-2}"
STORE_STEP="${write_interval:-1}"

PHI="${phi:-1e-2}"
BOX_AR="${aspect_ratio:-1}"
PCT_NONMAG="$( [ "$MIX_TYPES" = "TRUE" ] && echo "${percentnonmagpart:-0.5}" || echo 0 )"

# Parâmetros físicos
LAMBDA_VAL="${lambda:-1}"
ALPHA_VAL="${alpha_field:-0}"
PE_B="$( [ "$BROWNIAN" = "TRUE" ] && echo "${pecletnumber:-1}" || echo 0 )"
STK_T="$( [ "$PARTICLE_INERTIA" = "TRUE" ] && echo "${stokesnumber:-1e-1}" || echo 0 )"

# Flags padrão + as que você pediu menu
STATIC_SIM="$( [ "${choice1_05:-1}" = "0" ] && echo TRUE || echo FALSE )"  # Monte Carlo -> STATIC
RECORD_POS=$(bool_from_choice "${choice_rec_pos:-1}")
RECORD_VEL=$(bool_from_choice "${choice_rec_vel:-1}")
RECORD_DIP=$(bool_from_choice "${choice_rec_dip:-1}")
OVITO_FMT=$(bool_from_choice "${choice_ovito:-0}") # default OVITO TRUE

ORDERED_ARR=$(bool_from_choice "${choice_ordered_arr:-1}")
INIT_SPH_AGG=$(bool_from_choice "${choice_init_sph:-0}")
ORDERED_DIPOLES=$(bool_from_choice "${choice_ordered_dip:-1}")

# --- Sanitize booleans (1 linha) para evitar CR/LF acidentais -----------------
for v in SOLVE_TORQUE MAG_PART STATIC_SIM SEDIMENTATION FLUCT_MODE BROWNIAN HDI \
         CONTINUE_OLD PARTICLE_INERTIA PERIODIC_TORQUES PERIODIC_FORCES \
         RECORD_POS RECORD_VEL RECORD_DIP OVITO_FMT ORDERED_ARR MIX_TYPES \
         INIT_SPH_AGG ORDERED_DIPOLES APPLY_FIELD OSC_SHEAR TURN_ON_SHEAR DYN_SHEAR
do
  eval "$v=\$(oneline \"\${$v}\")"
done

# ===================== gera simconfig.dat =====================================
cat > simconfig.dat <<EOF
***********************************************
*                                             *
*  SIMMSUS:  SIMULATION CONFIGURATION FILE    *
*                                             *
***********************************************

- SIMULATION LOGICAL VARIABLES

 SOLVE TORQUE EQUATION..................: $SOLVE_TORQUE
 MAGNETIC PARTICLES.....................: $MAG_PART
 STATIC (T) OR DYNAMIC (F) SIMULATION...: $STATIC_SIM
 SEDIMENTATION..........................: $SEDIMENTATION
 FLUCTUATION MODE.......................: $FLUCT_MODE
 BROWNIAN MOTION........................: $BROWNIAN
 ACCOUNT HYDRODYNAMIC INTERACTIONS......: $HDI
 CONTINUE AN OLD SIMULATION.............: $CONTINUE_OLD
 PARTICLE INERTIA.......................: $PARTICLE_INERTIA

- DIPOLAR INTERACTION CALCULATIONS
 
 PERIODIC MAGNETIC TORQUES..............: $PERIODIC_TORQUES
 PERIODIC MAGNETIC FORCES...............: $PERIODIC_FORCES

- RECORD DATA OPTIONS
 
 RECORD POSITION IN FILE................: $RECORD_POS
 RECORD VELOCITY IN FILE................: $RECORD_VEL
 RECORD DIPOLE IN FILE..................: $RECORD_DIP
 OVITO(T) OR TECPLOT(F) OUTPUT FORMAT...: $OVITO_FMT

 - SIMULATION QUANTITIES

 NUMBER OF PARTICLES....................: $(zpad 5 "$N_PART")
 NUMBER OF REALIZATIONS.................: $(zpad 5 "$N_REALZ")

- INITIAL CONFIGURATION INFORMATION

 ORDERED (T) OR RANDOM (F) ARRANGEMENT..: $ORDERED_ARR
 MIX MAGNETIC AND NON MAGNETIC PARTICLES: $MIX_TYPES
 INITIAL SPHERICAL AGGREGATE............: $INIT_SPH_AGG
 ORDERED DIPOLES........................: $ORDERED_DIPOLES
 
- PARTICLE DISTRIBUTION DATA

 MONO (F) OR POLIDISPERSITY (T).........: $POLY
 VOLUME FRACTION OF PARTICLES...........: $(sci "$PHI")
 BOX ASPECT RATIO.......................: $(zpad 5 "$BOX_AR")
 PERCENTAGE OF NON-MAGNETIC PARTICLES...: $(sci "$PCT_NONMAG")

- MAGNETIC FIELD INFORMATION

 APPLY AN EXTERNAL MAGNETIC FIELD.......: $APPLY_FIELD
 POSITION OF THE EXTERNAL FIELD.........: $(zpad 5 "$FIELD_POS")
 OSCILLATORY FIELD......................: $OSC_FIELD
 ROTATING FIELD.........................: $ROT_FIELD
 NON-LINEAR DUFFING FIELD EXCITATION....: $DUFFING
 DOUBLE FREQUENCY FIELD EXCITATION......: $DOUBLE_FREQ
 DYNAMICAL INCREASE OF FIELD FREQUENCY..: $DYN_FREQ
 FREQUENCY 1 OF THE MAGNETIC FIELD......: $(sci "$omega1")
 FREQUENCY 2 OF THE MAGNETIC FIELD......: $(sci "$omega2")
 C1 PARAMETER FOR DUFFING EXCITATION....: $(sci "$c1")
 C2 PARAMETER FOR DUFFING EXCITATION....: $(sci "$c2")
 C3 PARAMETER FOR DUFFING EXCITATION....: $(sci "$c3")
 C4 PARAMETER FOR DUFFING EXCITATION....: $(sci "$c4")
 MAX FREQUENCY FOR DYNAMICAL INCREASE...: $(sci "$omega_max")
 NUMBER OF INTERVALS FOR DYN. INCREASE..: $(zpad 5 "$n_intervals")
 
- SHEAR RATE INFORMATION

 TURN ON SHEAR RATE.....................: $TURN_ON_SHEAR
 OSCILLATORY SHEAR......................: $OSC_SHEAR
 DYNAMICAL INCREASE OF SHEAR RATE.......: $DYN_SHEAR
 DIMENSIONLESS SHEAR RATE...............: $(sci "$GAMMA_D")
 FREQUENCY FOR THE OSCILLATORY SHEAR....: $(sci "$FREQ_SHEAR")

- PHYSICAL PARAMETERS

 LAMBDA ................................: $(sci "$LAMBDA_VAL")
 ALPHA .................................: $(sci "$ALPHA_VAL")
 BROWNIAN PECLET NUMBER.................: $(sci "$PE_B")
 TRANSLATIONAL STOKES NUMBER............: $(sci "$STK_T")

- NUMERICAL DATA

 SIMULATION TIME........................: $(sci "$T_SIM")
 NUMERICAL TIME-STEP....................: $(sci "$DT")
 STEP FOR STORING THE RESULTS...........: $(zpad 5 "$STORE_STEP")
 CONTINUE FROM ITERACTION NUMBER........: ${CONT_FROM:-10000}

- OPTIONAL DATA TREATMENT

 STATISTICAL ANALYSIS...................: FALSE
 CALCULATE THE STRUCTURE FACTOR.........: FALSE
 PRINT LOCAL MAPS OF PHI................: FALSE 
EOF

echo "✅ simconfig.dat gerado na pasta atual."
