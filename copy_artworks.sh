#!/usr/bin/env bash
rm -rf ./manuscript/artworks/primates/; rm -rf ./manuscript/artworks/mammals/; rm -rf ./manuscript/artworks/isopods/; rm -rf ./manuscript/artworks/drosophila/; rm -rf ./manuscript/artworks/simulations/
mkdir -p ./manuscript/artworks/primates/; mkdir ./manuscript/artworks/mammals/; mkdir ./manuscript/artworks/isopods/; mkdir ./manuscript/artworks/drosophila/; mkdir ./manuscript/artworks/simulations/

### Simulations ###
# Profiles correlation (png) #
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_profiles/np-5000-30_prefs.pdf ./manuscript/artworks/simulations/Profiles_SimuDiv.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_profiles/SimuDiv-node-False-1-read_prefs.pdf ./manuscript/artworks/simulations/Profiles_SimuDiv_Estimated.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_aa_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSel_ProfileCorrelation.png
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_aa_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMdutSel_ProfileCorrelation.png
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSelBranchNe_ProfileCorrelation.png
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_node_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMutSelBranchNe_ProfileCorrelation.png
# Branch correlation (pdf) #
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuFold_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuFold_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuFold_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuFold_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuFold_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuFold_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuFold_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuFold_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuFold_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuFold_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuDiv_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuPoly_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuGeo_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuGeo_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuGeo_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuGeo_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuGeo_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuGeo_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuGeo_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuGeo_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuGeo_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuGeo_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuFold_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuFold_BranchOmega_BranchCorrelation_LogOmega.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuDiv_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuDiv_BranchOmega_BranchCorrelation_LogOmega.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuPoly_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuPoly_BranchOmega_BranchCorrelation_LogOmega.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000BranchWise/inference_SimuGeo_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/BranchWise_SimuGeo_BranchOmega_BranchCorrelation_LogOmega.pdf
# ------------- #
# Profiles correlation (png) #
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_profiles/np-5000-498_prefs.pdf ./manuscript/artworks/simulations/Profiles_SimuPoly.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_profiles/SimuPoly-node-False-1-read_prefs.pdf ./manuscript/artworks/simulations/Profiles_SimuPoly_Estimated.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_aa_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/SimuDiv_SiteMutSel_ProfileCorrelation.png
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_aa_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/SimuPoly_SiteMutSel_ProfileCorrelation.png
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_node_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/SimuDiv_SiteMutSelBranchNe_ProfileCorrelation.png
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_profiles/correlation.aa-preferences-Simulation-1.png ./manuscript/artworks/simulations/SimuPoly_SiteMutSelBranchNe_ProfileCorrelation.png
# Branch correlation (pdf) #
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuFold_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/SimuFold_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuFold_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/SimuFold_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuFold_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/SimuFold_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuFold_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/SimuFold_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuFold_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/SimuFold_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/SimuDiv_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/SimuDiv_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/SimuDiv_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/SimuDiv_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/SimuDiv_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_trees/correlation.BranchTime-Simulation-1.pdf ./manuscript/artworks/simulations/SimuPoly_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_trees/correlation.LogPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/SimuPoly_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_trees/correlation.ContrastPopulationSize-Simulation-1.pdf ./manuscript/artworks/simulations/SimuPoly_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_trees/correlation.Log10BranchLength-Simulation-1.pdf ./manuscript/artworks/simulations/SimuPoly_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_node_trees/correlation.LogMutationRatePerTime-Simulation-1.pdf ./manuscript/artworks/simulations/SimuPoly_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuGeo_node_trees/correlation.BranchTime-Simulation-2.pdf ./manuscript/artworks/simulations/SimuGeo_SiteMutSelBranchNe_BranchCorrelation_BranchTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuGeo_node_trees/correlation.LogPopulationSize-Simulation-2.pdf ./manuscript/artworks/simulations/SimuGeo_SiteMutSelBranchNe_BranchCorrelation_LogPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuGeo_node_trees/correlation.ContrastPopulationSize-Simulation-2.pdf ./manuscript/artworks/simulations/SimuGeo_SiteMutSelBranchNe_BranchCorrelation_ContrastPopulationSize.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuGeo_node_trees/correlation.Log10BranchLength-Simulation-2.pdf ./manuscript/artworks/simulations/SimuGeo_SiteMutSelBranchNe_BranchCorrelation_Log10BranchLength.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuGeo_node_trees/correlation.LogMutationRatePerTime-Simulation-2.pdf ./manuscript/artworks/simulations/SimuGeo_SiteMutSelBranchNe_BranchCorrelation_LogMutationRatePerTime.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuFold_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/SimuFold_BranchOmega_BranchCorrelation_LogOmega.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuDiv_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/SimuDiv_BranchOmega_BranchCorrelation_LogOmega.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuPoly_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/SimuPoly_BranchOmega_BranchCorrelation_LogOmega.pdf
\cp -f ./DataSimulated/Experiments/Gen5Ma150L5000/inference_SimuGeo_nodeomega_trees/correlation.LogOmega-Simulation-1.pdf ./manuscript/artworks/simulations/SimuGeo_BranchOmega_BranchCorrelation_LogOmega.pdf

# Identifiability
\cp -f ./DataSimulated/Analysis/identifiability.Gen5Ma150L5000.pdf ./manuscript/artworks/simulations/identifiability.pdf

### Isopods ###
# Experiment correlation (pdf) #
\cp -f ./DataEmpirical/Analysis/identifiability.Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6.pdf ./manuscript/artworks/isopods/identifiability.pdf

\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize-1-2.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-2_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize-1-3.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-3_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize-1-4.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-4_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize-1-5.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-5_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize-1-6.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-6_LogPopulationSize.pdf

\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize-1-2.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-2_ContrastPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize-1-3.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-3_ContrastPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize-1-4.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-4_ContrastPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize-1-5.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-5_ContrastPopulationSize.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize-1-6.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-6_ContrastPopulationSize.pdf

\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogMutationRatePerTime-1-2.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-2_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogMutationRatePerTime-1-3.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-3_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogMutationRatePerTime-1-4.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-4_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogMutationRatePerTime-1-5.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-5_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogMutationRatePerTime-1-6.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-6_LogMutationRatePerTime.pdf

\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.Log10BranchLength-1-2.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-2_Log10BranchLength.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.Log10BranchLength-1-3.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-3_Log10BranchLength.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.Log10BranchLength-1-4.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-4_Log10BranchLength.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.Log10BranchLength-1-5.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-5_Log10BranchLength.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.Log10BranchLength-1-6.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-6_Log10BranchLength.pdf

\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.BranchTime-1-2.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-2_BranchTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.BranchTime-1-3.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-3_BranchTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.BranchTime-1-4.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-4_BranchTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.BranchTime-1-5.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-5_BranchTime.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.BranchTime-1-6.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep-1-6_BranchTime.pdf
# Trees (pdf) #
\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id0/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R1_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id1/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R2_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id2/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R3_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id3/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R4_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id4/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R5_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id5/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R6_LogPopulationSize.pdf

\cp -f ./DataEmpirical/Experiments/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_Id0/inference_node_trees/node_False_1_run.LogMutationRatePerTime.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_R1_LogMutationRatePerTime.pdf
# Stats (pdf) #
#\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize.nhx.eco.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize_eco.pdf
#\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize.nhx.eye.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize_eye.pdf
#\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize.nhx.pig.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize_pig.pdf
#\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize.nhx.eco.merged.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize_eco_merged.pdf
#\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize.nhx.eye.merged.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize_eye_merged.pdf
#\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.ContrastPopulationSize.nhx.pig.merged.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize_pig_merged.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize.nhx.eco.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize_eco.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize.nhx.eye.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize_eye.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize.nhx.pig.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize_pig.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize.nhx.eco.merged.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize_eco_merged.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize.nhx.eye.merged.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize_eye_merged.pdf
\cp -f ./DataEmpirical/Analysis/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6/node_False_1_run.LogPopulationSize.nhx.pig.merged.pdf ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize_pig_merged.pdf
\cp -f ./DataEmpirical/Analysis/DataFrame/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_node_False_1_run.LogPopulationSize.tex ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogPopulationSize.tex
\cp -f ./DataEmpirical/Analysis/DataFrame/Ncat50_Isopods_rootedtree.nhx_cds.highcoverage.list_Sample12_Replicates6_node_False_1_run.LogMutationRatePerTime.tex ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_Rep_LogMutationRatePerTime.tex

\cp -f ./DataEmpirical/Analysis/DataFrame/anova_ocular_structure.txt ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_anova_ocular_structure.txt
\cp -f ./DataEmpirical/Analysis/DataFrame/anova_pigmentation.txt ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_anova_pigmentation.txt
\cp -f ./DataEmpirical/Analysis/DataFrame/anova_habitat.txt ./manuscript/artworks/isopods/12CDS_SiteMutSelBranchNe_anova_habitat.txt
### Mammals ###
# Experiment correlation (pdf) #
\cp -f ./DataEmpirical/Analysis/identifiability.Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4.pdf ./manuscript/artworks/mammals/identifiability.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.LogPopulationSize-1-2.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogPopulationSize-1-2.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.LogPopulationSize-1-3.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogPopulationSize-1-3.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.LogPopulationSize-1-4.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogPopulationSize-1-4.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.LogMutationRatePerTime-1-2.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogMutationRatePerTime-1-2.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.LogMutationRatePerTime-1-3.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogMutationRatePerTime-1-3.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.LogMutationRatePerTime-1-4.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogMutationRatePerTime-1-4.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.ContrastPopulationSize-1-2.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize-1-2.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.ContrastPopulationSize-1-3.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize-1-3.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.ContrastPopulationSize-1-4.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_ContrastPopulationSize-1-4.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.Log10BranchLength-1-2.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_Log10BranchLength-1-2.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.Log10BranchLength-1-3.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_Log10BranchLength-1-3.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.Log10BranchLength-1-4.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_Log10BranchLength-1-4.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.BranchTime-1-2.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_BranchTime-1-2.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.BranchTime-1-3.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_BranchTime-1-3.pdf
\cp -f ./DataEmpirical/Analysis/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4/node_False_1_run.BranchTime-1-4.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_BranchTime-1-4.pdf
\cp -f ./DataEmpirical/Analysis/DataFrame/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_node_False_1_run.LogPopulationSize.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogPopulationSize.tex
\cp -f ./DataEmpirical/Analysis/DataFrame/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_node_False_1_run.LogMutationRatePerTime.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_Rep_LogMutationRatePerTime.tex

# Profiles correlation (png) #
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_aa_profiles/correlation.aa-preferences-1-2.png ./manuscript/artworks/mammals/18CDS_SiteMutSel_R1_ProfileCorrelation.png
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_profiles/correlation.aa-preferences-1-2.png ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_ProfileCorrelation.png
# Trees (pdf) #
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id1/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R2_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id2/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R3_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id3/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R4_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_trees/node_False_1_run.LogMutationRatePerTime.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_trees/node_False_1_run.TraitsAdult_weight_.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_LogAdult_weight.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_trees/node_False_1_run.TraitsFemale_maturity_.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_LogFemale_maturity.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_trees/node_False_1_run.TraitsMaximum_longevity_.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_LogMaximum_longevity.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_nodeomega_trees/nodeomega_False_1_run.LogMutationRatePerTime.pdf ./manuscript/artworks/mammals/18CDS_BranchOmega_R1_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_nodeomega_trees/nodeomega_False_1_run.LogOmega.pdf ./manuscript/artworks/mammals/18CDS_BranchOmega_R1_LogdNdS.pdf
# Traits correlation matrix (pdf & tex) #
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/inference_node_trees/correlation.LogPopulationSize-1-2.pdf ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_LogPopulationSizeCorrelation.pdf
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/CorrelationMatrices/node_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/CorrelationMatrices/node_False_1_covariance.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/CorrelationMatrices/node_False_1_partial_coefficient.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R1_TraitsPartialCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id1/CorrelationMatrices/node_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R2_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id2/CorrelationMatrices/node_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R3_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id3/CorrelationMatrices/node_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_SiteMutSelBranchNe_R4_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/CorrelationMatrices/nodeomega_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R1_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/CorrelationMatrices/nodeomega_False_1_covariance.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R1_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id0/CorrelationMatrices/nodeomega_False_1_partial_coefficient.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R1_TraitsPartialCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id1/CorrelationMatrices/nodeomega_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R2_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id1/CorrelationMatrices/nodeomega_False_1_covariance.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R2_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id1/CorrelationMatrices/nodeomega_False_1_partial_coefficient.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R2_TraitsPartialCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id2/CorrelationMatrices/nodeomega_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R3_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id2/CorrelationMatrices/nodeomega_False_1_covariance.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R3_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id2/CorrelationMatrices/nodeomega_False_1_partial_coefficient.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R3_TraitsPartialCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id3/CorrelationMatrices/nodeomega_False_1_correlation.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R4_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id3/CorrelationMatrices/nodeomega_False_1_covariance.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R4_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/Cat50_OrthoMam_rootedtree.lht.nhx_cds.highcoverage.list_Sample18_Replicates4_Id3/CorrelationMatrices/nodeomega_False_1_partial_coefficient.tex ./manuscript/artworks/mammals/18CDS_BranchOmega_R4_TraitsPartialCorrelation.tex
### Primates ###
# Profiles correlation (png) #
\cp -f ./DataEmpirical/Experiments/IWCat50LHT_Primates_CDS.ali_rootedtree.nhx/inference_aa_profiles/correlation.aa-preferences-1-2.png ./manuscript/artworks/primates/SiteMutSel_ProfileCorrelation.png
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_profiles/correlation.aa-preferences-1-2.png  ./manuscript/artworks/primates/SiteMutSelBranchNe_ProfileCorrelation.png
# Trees (pdf)
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/correlation.LogPopulationSize-1-2.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_LogPopulationSizeCorrelation.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.LogMutationRatePerTime.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.LogPopulationSize.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_LogPopulationSize.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.Traitsgeneration_time.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_Loggeneration_time.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.Traitslongevity.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_Loglongevity.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.Traitsmass.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_Logmass.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.Traitsmaturity.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_Logmaturity.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.TraitspiNpiS.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_LogpiNpiS.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_node_trees/node_False_1_run.TraitspiS.pdf ./manuscript/artworks/primates/SiteMutSelBranchNe_LogpiS.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_nodeomega_trees/nodeomega_False_1_run.LogMutationRatePerTime.pdf ./manuscript/artworks/primates/BranchOmega_LogMutationRatePerTime.pdf
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/inference_nodeomega_trees/nodeomega_False_1_run.LogOmega.pdf ./manuscript/artworks/primates/BranchOmega_LogOmega.pdf
# Traits correlation matrix (pdf & tex)
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/CorrelationMatrices/node_False_1_correlation.tex ./manuscript/artworks/primates/SiteMutSelBranchNe_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/CorrelationMatrices/node_False_1_covariance.tex ./manuscript/artworks/primates/SiteMutSelBranchNe_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/CorrelationMatrices/node_False_1_partial_coefficient.tex ./manuscript/artworks/primates/SiteMutSelBranchNe_TraitsPartialCorrelation.tex
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/CorrelationMatrices/nodeomega_False_1_correlation.tex ./manuscript/artworks/primates/BranchOmega_TraitsCorrelation.tex
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/CorrelationMatrices/nodeomega_False_1_covariance.tex ./manuscript/artworks/primates/BranchOmega_TraitsCovariance.tex
\cp -f ./DataEmpirical/Experiments/IWCat50LHTCalib_Primates_CDS.ali_rootedtree.nhx/CorrelationMatrices/nodeomega_False_1_partial_coefficient.tex ./manuscript/artworks/primates/BranchOmega_TraitsPartialCorrelation.tex
