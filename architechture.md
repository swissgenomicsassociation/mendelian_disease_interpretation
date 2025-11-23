# Architecture plan for the consensus guideline

## **0. Purpose and scope**

This guideline defines the national minimum contract for evidence based variant interpretation in Mendelian disease. Its central output is a structured set of machine readable evidence flags and the rules that generate them. These flags represent falsification tests applied to every variant and provide a consistent evidence state that can be compared, exchanged, and audited across institutions, platforms, and time.

The framework follows reverse reasoning. For every variant proposed by an upstream tool, the system evaluates each evidence domain for signals that contradict, weaken, or fail to support the hypothesis that the variant explains disease. All evidence is recorded as explicit flags rather than implicit assumptions.

The guideline describes how evidence is represented and interpreted and does not prescribe sequencing workflows, platforms, pipelines, or prioritisation tools. Laboratories and companies may use any methods they choose, provided they supply the minimal identifiers and evidence required for interoperability.

## **1. Reframing the guideline: from manual protocols to computable standards**

Traditional genetics guidelines describe manual steps that experts are expected to follow. These documents list actions, recommendations, and interpretive tasks that rely on specialist judgement and remain vulnerable to variation across laboratories. Such documents cannot ensure reproducibility at national scale, and they cannot be executed by automated systems or validated computationally.

Our approach removes this dependency. Instead of prescribing manual protocols, the guideline defines a set of **computable evidence rules** and the **flags** that these rules generate. Each rule tests whether a specific domain of evidence supports, weakens, or contradicts the hypothesis that a variant is the cause of disease. Every rule is implemented in a structured, machine readable format and is versioned, auditable, and openly released.

This creates several advantages:

* **Rules replace prose**. Interpretation steps are encoded as formal logic rather than narrative instruction.
* **Flags replace subjective judgement**. Each domain outputs an explicit, binary or categorical evidence state that can be reproduced exactly by any institution.
* **Public versioning replaces expert drift**. Rule sets evolve through transparent releases. The most widely adopted versions naturally become the stable national standard.
* **Automation replaces manual review**. Upstream tools may vary, but all downstream evaluation is performed through the same falsification framework.
* **Non expert misuse is prevented**. Laboratories and companies cannot bypass or reinterpret the logic because the system defines what must be computed, not how a human should behave.
* **National harmonisation emerges organically**. Shared rule sets allow institutions to align without imposing specific pipelines, vendors, or technologies.

This strategy ensures that the national evidence model remains computable, interoperable, and future proof. It aligns with digital health architecture principles, avoids reliance on scarce manual expertise, and transforms variant interpretation from a craft into a reproducible scientific procedure grounded in explicit rules and counter evidence.

## **2. Architectural context**

National genomic infrastructure must be durable while remaining adaptable to scientific change. Clinical records, identifiers, and audit trails require stable structure, consistent schemas, and strong governance. Biological meaning, phenotype models, and evidence resources evolve continuously and must be represented in flexible, semantic formats.

The guideline is designed for this dual environment. Stable components define what must be recorded to ensure traceability and reproducibility. Flexible components define how evidence rules and knowledge sources can be updated without invalidating past interpretations. Clear separation between these layers ensures that future reference data, ontologies, or annotation models can be incorporated without disrupting national interoperability.

## **3. Position within the three pillar framework**

The Swiss Genomics Association maintains a national architecture consisting of three coordinated pillars. Pillar 1 defines sample and sequencing provenance. Pillar 2 provides normalised analysis variables and identifiers. Pillar 3 defines computable evidence and interpretation.

The present work forms the foundation of Pillar 3. It specifies the falsification logic, evidence rules, and flag outputs that operate on the harmonised inputs defined in Pillar 2 and the provenance captured in Pillar 1. The guideline does not reproduce the full three pillar model but relies on it for complete interoperability.

## **4. System design principles**

The guideline is built on principles that support national scale genomics and long term interpretive consistency.

* Workflow independence. Any institution can adopt the standard without modifying local pipelines.
* Machine readable evidence. All rules are encoded, versioned, and auditable.
* Explicit uncertainty. Absence of evidence is recorded directly and never inferred.
* Full lineage. Every flag is traceable to the rule, reference resource, and input data that produced it.
* Decentralised use. Laboratories and companies maintain their own systems while sharing a common evidence contract.
* Alignment with SPHN, GA4GH, RDF semantics, and internationally recognised data standards.
* Versioned schemas, public reference implementations, and forward compatible rule updates.
* Clinical traceability. Each interpretation must be reconstructable from the evidence state.

## **5. Implemented rule set**

The current reference implementation provides a reproducible set of functions driven by a central YAML rule file. These functions operate independently of any variant calling or annotation pipeline.

* **01_qv_logic.R**
  Applies the complete evidence rule set encoded in YAML.

* **00_MAX_FREQ_function.R**
  Computes population frequency flags for common, rare, ultrarare, and missing states using gnomAD and related references.

* **00_import_exomiser.R**
  Imports Exomiser variant files, normalises identifiers, and constructs the base variant table.

* **00_MOI_function.R**
  Evaluates segregation and inheritance compatibility. Derives trio genotypes, standardises REF, HET, HOM, and HEMI, and assesses AD, AR, and XL disease mechanisms. Flags consistent, conflicting, and incomplete states.

* **00_UNIPROT_function.R**
  Maps UniProt features to genomic positions, detects domain overlap, truncation, and predicted NMD, and records functional evidence states.

* **99_generate_interpretations.R**
  Experimental LLM based summarisation. Not part of the normative guideline.

## **6. Rules under development**

Additional domains will expand the falsification framework. Each rule contributes explicit evidence flags and does not remove variants.

* OMIM disease consistency
* ClinGen gene validity
* HPO phenotype alignment
* Transcript concordance for MANE
* Splicing proximity and prediction evidence
* Read level quality and artefact risk
* Presence in unrelated controls
* ClinVar conflict
* HGMD support
* Literature concordance
* Loss of function context using LOEUF and LOFTee
* Disease spectrum fit for severity and organ system involvement

## **7. Evidence domains**

Each evidence domain tests whether available data contradict the causal hypothesis.

### **7.1 Functional and molecular evidence**

UniProt features are mapped to variant coordinates at amino acid resolution. Rules test whether a variant overlaps a functional region, removes downstream domains, or triggers predicted NMD. All outcomes are encoded as evidence flags.

### **7.2 Population frequency**

Population based rules evaluate whether the allele frequency supports or contradicts a Mendelian disease hypothesis. Flags represent rarity, ultrararity, commonness, or missingness.

### **7.3 Mode of inheritance**

Segregation checks test whether observed genotypes align with the expected disease mechanism. Missing parental data, incompatible patterns, and unexpected zygosity states are flagged explicitly. All logic relies on standardised genotype encodings.

## **8. Final flag set**

The current complete flag set is:

```
flag_gt_valid
flag_moi_parent_gt_missing_mother
flag_moi_parent_gt_missing_father
flag_moi_parent_gt_missing_any
flag_moi_parent_gt_hom_mother
flag_moi_parent_gt_hom_father
flag_moi_parent_gt_hom_any
flag_moi_parent_conflict_AD
flag_moi_parent_conflict_AR
flag_moi_parent_conflict_XR
flag_moi_parent_conflict_any
flag_popfreq_common
flag_popfreq_rare
flag_popfreq_ultrarare
flag_missing_popfreq
flag_uniprot_hits_any_feature
flag_uniprot_hits_domain_like
flag_uniprot_hits_structural_like
flag_uniprot_hits_ptm_like
flag_uniprot_hits_binding_like
flag_uniprot_hits_variant_like
flag_uniprot_is_lof
flag_uniprot_predicted_nmd
flag_uniprot_truncates_feature
```

These flags define the reproducible evidence state for each variant and form the core input to Pillar 3 probabilistic interpretation.


Here is the conceptual reformulation you are reaching for, written with clarity, precision, and in full alignment with your digital, rule based national architecture.

This replaces the outdated “list of things to do” model seen in older genetics guidelines with a modern specification driven by machine readable rules, versioned evidence states, and reproducible computation.

You can drop this directly into your architecture overview, methods section, or as a framing paragraph for the full guideline.


