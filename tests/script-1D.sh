#!/bin/bash

# 1D tests
rm -rf 1D/*/output-tabs-SE27XX/ 
rm -rf 1D/*/geotop.log-SE27XX

# (1) B2_BeG_017
../../geotop_2.0/bin/geotop-2.0.0 1D/B2_BeG_017 1>>1D/B2_BeG_017/stdout.SE27XX 2>>1D/B2_BeG_017/stderr.SE27XX
mv 1D/B2_BeG_017/geotop.log 1D/B2_BeG_017/geotop.log-SE27XX
mv 1D/B2_BeG_017/output-tabs 1D/B2_BeG_017/output-tabs-SE27XX
mkdir 1D/B2_BeG_017/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (2) Bro
../../geotop_2.0/bin/geotop-2.0.0 1D/Bro 1>>1D/Bro/stdout.SE27XX 2>>1D/Bro/stderr.SE27XX
mv 1D/Bro/geotop.log 1D/Bro/geotop.log-SE27XX
mv 1D/Bro/output-tabs 1D/Bro/output-tabs-SE27XX
mkdir 1D/Bro/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (3) Calabria
../../geotop_2.0/bin/geotop-2.0.0 1D/Calabria 1>>1D/Calabria/stdout.SE27XX 2>>1D/Calabria/stderr.SE27XX
mv 1D/Calabria/geotop.log 1D/Calabria/geotop.log-SE27XX
mv 1D/Calabria/output-tabs 1D/Calabria/output-tabs-SE27XX
mkdir 1D/Calabria/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (4) ColdelaPorte
../../geotop_2.0/bin/geotop-2.0.0 1D/ColdelaPorte 1>>1D/ColdelaPorte/stdout.SE27XX 2>>1D/ColdelaPorte/stderr.SE27XX
mv 1D/ColdelaPorte/geotop.log 1D/ColdelaPorte/geotop.log-SE27XX
mv 1D/ColdelaPorte/output-tabs 1D/ColdelaPorte/output-tabs-SE27XX
mkdir 1D/ColdelaPorte/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (5) CostantMeteo
../../geotop_2.0/bin/geotop-2.0.0 1D/CostantMeteo 1>>1D/CostantMeteo/stdout.SE27XX 2>>1D/CostantMeteo/stderr.SE27XX
mv 1D/CostantMeteo/geotop.log 1D/CostantMeteo/geotop.log-SE27XX
mv 1D/CostantMeteo/output-tabs 1D/CostantMeteo/output-tabs-SE27XX
mkdir 1D/CostantMeteo/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (6) Jungfraujoch
../../geotop_2.0/bin/geotop-2.0.0 1D/Jungfraujoch 1>>1D/Jungfraujoch/stdout.SE27XX 2>>1D/Jungfraujoch/stderr.SE27XX
mv 1D/Jungfraujoch/geotop.log 1D/Jungfraujoch/geotop.log-SE27XX
mv 1D/Jungfraujoch/output-tabs 1D/Jungfraujoch/output-tabs-SE27XX
mkdir 1D/Jungfraujoch/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (7) Matsch_B2_Ref_007
../../geotop_2.0/bin/geotop-2.0.0 1D/Matsch_B2_Ref_007 1>>1D/Matsch_B2_Ref_007/stdout.SE27XX 2>>1D/Matsch_B2_Ref_007/stderr.SE27XX
mv 1D/Matsch_B2_Ref_007/geotop.log 1D/Matsch_B2_Ref_007/geotop.log-SE27XX
mv 1D/Matsch_B2_Ref_007/output-tabs 1D/Matsch_B2_Ref_007/output-tabs-SE27XX
mkdir 1D/Matsch_B2_Ref_007/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (8) Matsch_P2_Ref_007
../../geotop_2.0/bin/geotop-2.0.0 1D/Matsch_P2_Ref_007 1>>1D/Matsch_P2_Ref_007/stdout.SE27XX 2>>1D/Matsch_P2_Ref_007/stderr.SE27XX
mv 1D/Matsch_P2_Ref_007/geotop.log 1D/Matsch_P2_Ref_007/geotop.log-SE27XX
mv 1D/Matsch_P2_Ref_007/output-tabs 1D/Matsch_P2_Ref_007/output-tabs-SE27XX
mkdir 1D/Matsch_P2_Ref_007/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (9) PureDrainage
../../geotop_2.0/bin/geotop-2.0.0 1D/PureDrainage 1>>1D/PureDrainage/stdout.SE27XX 2>>1D/PureDrainage/stderr.SE27XX
mv 1D/PureDrainage/geotop.log 1D/PureDrainage/geotop.log-SE27XX
mv 1D/PureDrainage/output-tabs 1D/PureDrainage/output-tabs-SE27XX
mkdir 1D/PureDrainage/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (10) PureDrainageFaked
../../geotop_2.0/bin/geotop-2.0.0 1D/PureDrainageFaked 1>>1D/PureDrainageFaked/stdout.SE27XX 2>>1D/PureDrainageFaked/stderr.SE27XX
mv 1D/PureDrainageFaked/geotop.log 1D/PureDrainageFaked/geotop.log-SE27XX
mv 1D/PureDrainageFaked/output-tabs 1D/PureDrainageFaked/output-tabs-SE27XX
mkdir 1D/PureDrainageFaked/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (11) PureDrainageRainy
../../geotop_2.0/bin/geotop-2.0.0 1D/PureDrainageRainy 1>>1D/PureDrainageRainy/stdout.SE27XX 2>>1D/PureDrainageRainy/stderr.SE27XX
mv 1D/PureDrainageRainy/geotop.log 1D/PureDrainageRainy/geotop.log-SE27XX
mv 1D/PureDrainageRainy/output-tabs 1D/PureDrainageRainy/output-tabs-SE27XX
mkdir 1D/PureDrainageRainy/output-tabs
# ----------------------------------------------------------------------------------------------------------------
# (12) PureDrainageRainySlope
../../geotop_2.0/bin/geotop-2.0.0 1D/PureDrainageRainySlope 1>>1D/PureDrainageRainySlope/stdout.SE27XX 2>>1D/PureDrainageRainySlope/stderr.SE27XX
mv 1D/PureDrainageRainySlope/geotop.log 1D/PureDrainageRainySlope/geotop.log-SE27XX
mv 1D/PureDrainageRainySlope/output-tabs 1D/PureDrainageRainySlope/output-tabs-SE27XX
mkdir 1D/PureDrainageRainySlope/output-tabs
# ----------------------------------------------------------------------------------------------------------------
echo  1D/*/output-tabs/ | xargs -n 1 cp .placeholder
echo  1D/*/output-tabs-SE27XX/ | xargs -n 1 cp .placeholder

# ls -l 1D/ | awk '{print $9}'
