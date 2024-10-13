# chmod +x run_noisy_pixels.sh
#!/bin/bash

# A few useful timestamps [CET]:
# 1714514400 = 01/05/2024, 00:00:00
# 1719784800 = 01/07/2024, 00:00:00 
# 1725141600 = 01/09/2024, 00:00,00
# 1730415600 = 01/11/2024, 00:00:00 

# UNIX timestamp (type long, second precision -> 10 digits):
trg_start_min=1714514400
trg_start_max=1730415600

#root -q 'noisy_pixels_count.cxx+('$trg_start_min','$trg_start_max')'
root -q 'noisy_pixels_plots.cxx('$trg_start_min','$trg_start_max')'