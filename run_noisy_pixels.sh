# chmod +x run_noisy_pixels.sh
#!/bin/bash

# A few useful timestamps:
# 1717192800 = 01/06/2024, 00:00:00 
# 1723672800 = 15/08/2024, 00:00,00
# 1728424800 = 09/10/2024, 00:00:00 

# UNIX timestamp (type long, second precision -> 10 digits):
trg_start_min=1717192800 
trg_start_max=1723672800

root -q 'noisy_pixels.cxx('$trg_start_min','$trg_start_max')'