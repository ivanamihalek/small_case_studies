# MYO7A gene therapy veriants

The **MYO7A variants most amenable to gene augmentation therapy are biallelic loss-of-function (LoF) variants**, particularly nonsense, frameshift and splice-disrupting variants that result in absent or severely reduced functional myosin VIIA.

However, there is an important distinction: **MYO7A gene therapy is largely mutation-independent.** Unlike exon-skipping or antisense approaches, gene augmentation introduces a functional copy of the gene rather than correcting the patient's specific mutation.

Consequently, the most important determinants of therapeutic amenability are the molecular mechanism, residual protein function, target-cell viability and ability to deliver the therapeutic construct, rather than the precise variant position.

## 1. Ranking MYO7A variant classes by gene augmentation suitability

| Variant class | Predicted amenability | Rationale |
|---|---|---|
| Nonsense / premature stop | Very high | Often results in NMD or truncated, non-functional protein. |
| Frameshift | Very high | Usually causes loss of functional MYO7A. |
| Canonical splice-site variants | Very high | Frequently result in aberrant splicing and LoF. |
| Large deletions | Very high | Particularly suitable for replacement if they eliminate functional protein. |
| Missense variants causing protein instability | High | Wild-type MYO7A could restore functional protein levels. |
| Missense variants affecting motor or functional domains | High, with caveats | Replacement may work, but residual or abnormal protein function needs consideration. |
| In-frame deletions / insertions | Moderate to high | Depends on the functional consequences of the altered protein. |
| Variants with substantial residual function | Potentially high | Additional functional MYO7A may improve cellular function, but this requires experimental validation. |
| Dominant-negative or gain-of-function variants | Lower / uncertain | Adding wild-type protein might not overcome the pathogenic effect of the mutant protein. |

One caveat: this is a mechanistic ranking, not a clinically validated ranking of individual MYO7A variants.

## 2. Why variant position is less important than molecular mechanism

MYO7A encodes myosin VIIA, a large protein with several functional domains:

- N-terminal motor domain: ATPase activity and actin interaction.
- Neck region: IQ motifs and calmodulin binding.
- Coiled-coil region: dimerisation.
- C-terminal tail: cargo binding and protein interactions, including two MyTH4 domains and FERM domains.

A missense variant in the motor domain might severely impair ATPase activity, whereas another missense variant in the same domain might retain substantial function.

For gene augmentation, this distinction is less important than it would be for a mutation-specific treatment.

Consider three scenarios:

1. **Homozygous nonsense variant:** The endogenous protein is absent or severely truncated. Supplying wild-type MYO7A is conceptually straightforward.
2. **Compound heterozygous missense variants:** Both alleles encode dysfunctional proteins. Wild-type MYO7A could theoretically restore function.
3. **Missense variant producing a dominant-negative protein:** The mutant protein may interfere with normal cellular processes. Adding wild-type MYO7A might not be sufficient.

Importantly, MYO7A-associated Usher syndrome is predominantly autosomal recessive. Dominant-negative mechanisms should therefore be investigated rather than assumed for individual variants.

A particularly relevant observation is that a study of 62 patients with biallelic MYO7A variants found no significant relationship between variant type or protein-domain location and visual disease severity or progression. This argues against selecting gene augmentation candidates based simply on variant position. [PubMed Central (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9540557/?utm_source=chatgpt.com)

## 3. The actual gene therapy challenge: MYO7A is too large for conventional AAV

This is arguably more important than the variant classification.

MYO7A cDNA is approximately 6.7 kb, exceeding the approximately 4.7 kb packaging capacity of conventional AAV vectors.

Several approaches have been investigated:

| Approach | Technology | Status / evidence |
|---|---|---|
| UshStat | EIAV lentiviral vector | Previously investigated in clinical trials |
| Dual AAV8-MYO7A | Split MYO7A cDNA, reconstituted intracellularly | Preclinical retinal efficacy |
| Dual AAV approaches | Split-gene delivery | Additional preclinical development |

The dual-AAV approach is particularly interesting because it allows delivery of the complete functional MYO7A coding sequence.

In 2023, Ferla et al. demonstrated improvement of retinal abnormalities in a mouse model using dual AAV8-MYO7A. They also evaluated vector safety and biodistribution in non-human primates. [PubMed Central (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9996380/?utm_source=chatgpt.com)

There is also evidence of improved vestibular function following dual-AAV MYO7A delivery in a shaker-1 mouse model. However, auditory function was not restored in that experiment. [PubMed Central (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10491803/?utm_source=chatgpt.com)

## 4. Clinical translation: an important development in 2026

There is now a particularly relevant clinical programme:

**AAVB-081 (AAVantgarde Bio)**

- Dual AAV8.MYO7A gene replacement.
- Subretinal administration.
- Phase I/II LUCE-1 clinical trial (NCT06591793).
- Designed for patients with MYO7A-associated Usher syndrome type 1B.

An important update: on 5 October 2026, AAVantgarde announced that all 15 participants had been treated. Preliminary findings included early signals of improved visual function in some participants, alongside encouraging initial safety observations. These are preliminary company-reported results, not evidence of established clinical efficacy. [AAVantgarde](https://www.aavantgarde.com/en/news/aavantgarde-presents-encouraging-clinical-data-from-aavb-081-and-preclinical-data-supporting-aavb-039-at-euretina-2026/?utm_source=chatgpt.com)

The trial is particularly interesting because its inclusion criteria are based on a molecular diagnosis of MYO7A-associated USH1B, rather than restricting eligibility to a particular mutation class.

## 5. What I would prioritise in a MYO7A variant dataset

If the objective is to identify patients most likely to benefit from gene augmentation, I would prioritise the following features:

| Feature | Priority | Interpretation |
|---|---|---|
| Biallelic pathogenic variants | Essential | Establishes the recessive disease mechanism. |
| Predicted LoF | High | Supports replacement rather than mutation-specific correction. |
| Nonsense-mediated decay | High | Supports absence of functional endogenous protein. |
| Absence of dominant-negative effects | High | Reduces concern about interference from mutant protein. |
| Retained photoreceptor structure | Very high | Indicates potentially treatable retinal cells remain. |
| Retained RPE integrity | High | Relevant to MYO7A-dependent retinal functions. |
| Variant location | Secondary | Less important for full-length gene replacement. |
| Residual MYO7A protein function | Important | Helps distinguish null alleles from hypomorphic alleles. |

One additional consideration is MYO7A alternative splicing. There are two major retinal isoforms, IF1 and IF2, with potentially overlapping but also distinct functions. Therapeutic construct design should account for both. [PubMed Central (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10984346/?utm_source=chatgpt.com)

---

## 6. My overall conclusion

For a gene augmentation programme, I would use this prioritisation:

1. **Highest priority:** Patients with biallelic nonsense, frameshift or severe splice-disrupting variants.
2. **High priority:** Patients with biallelic missense variants causing severe loss of protein function.
3. **Intermediate priority:** Patients with hypomorphic missense or in-frame variants.
4. **Lower or uncertain priority:** Patients with evidence of dominant-negative or gain-of-function mechanisms.

Crucially, I would not exclude missense variants simply because they retain some protein expression.

**The most interesting research question is whether residual MYO7A function predicts the magnitude of rescue following gene augmentation.** That is a more useful hypothesis to investigate than whether variants in one particular protein domain are more amenable than variants in another.

### Key papers and resources

1. [Genotype–phenotype correlation in patients with MYO7A variants (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9540557/)
2. [Dual AAV MYO7A gene therapy: retinal efficacy and safety (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9996380/)
3. [Dual AAV MYO7A gene therapy: vestibular rescue (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10491803/)
4. [MYO7A isoforms and implications for gene therapy (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10984346/)
5. [LUCE-1 clinical trial, NCT06591793](https://clinicaltrials.gov/study/NCT06591793)

**My recommendation:** If you're developing a computational prioritisation framework, distinguish *molecular amenability to gene augmentation* from *predicted clinical benefit*. They are related, but not interchangeable. A null variant can be an excellent molecular target while the patient has too little surviving photoreceptor tissue for meaningful visual rescue.
