# Coupling Europe's Capacity Markets
Julia implementation of the long-run equilibrium model described in:

Adekola, de Vries, Bruninx (2026). Coupling Europe's Capacity Markets. arXiv:2603.08248
https://arxiv.org/abs/2603.08248

## Abstract
European Member States are increasingly introducing national capacity mechanisms (CMs) to manage growing adequacy risks. However, isolated national CMs are inefficient in highly interconnected electricity systems, such as the European system. While progress has been made in facilitating cross-border participation by generation capacity in CMs, existing arrangements are prone to under- or over-investment and do not properly value the contribution of interconnection capacity to Member States' adequacy targets. In this paper, we propose a novel conceptual design for a coupled European capacity market that utilises the logic of flow-based market coupling. In a comparative analysis of different market design scenarios in an illustrative multi-zone case study, using a bespoke long-run equilibrium problem, we show that the proposed flow-based coupling of capacity markets reduces system costs by harnessing available capacity in neighbouring market zones while ensuring deliverability with respect to network constraints in all scarcity situations.

## Model
The model is a Mixed Complementarity Problem (MCP) solved via the Alternating Direction Method of Multipliers (ADMM). It includes:
Generation Agents: Investment and dispatch optimization.
Consumer Agents: Elastic and inelastic demand response.
Energy & Capacity Market Operators: Flow-based clearing and cross-border obligation management.


## Requirements
Julia - https://julialang.org/
Gurobi (with valid licence) - https://www.gurobi.com/

## Reference
If you use this code in your research, please cite:
Adekola, K., de Vries, L., & Bruninx, K. (2026). Coupling Europe's Capacity Markets. arXiv:2603.08248.

Disclaimer: The referenced paper and this code are currently under peer-review and may be subject to changes.
