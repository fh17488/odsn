#!/bin/bash

python3 ./collate_heatmap_data_YAvg.py small-world opinion_dynamics_SW_K128 K128 &
python3 ./collate_heatmap_data_YAvg.py random opinion_dynamics_ER_K128 K128 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_RL_K128 K128 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_RL_K16_D0 K16_D0 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_RL_K16_D0_25 K16_D0_25 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_K16_D0_5 K16_D0_5 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_K16_D0_75 K16_D0_75 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_K16_R0 K16_R0 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_K16_R0_25 K16_R0_25 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_K16_R0_75 K16_R0_75 &
python3 ./collate_heatmap_data_YAvg.py ring-lattice opinion_dynamics_K16_R1 K16_R1 &
