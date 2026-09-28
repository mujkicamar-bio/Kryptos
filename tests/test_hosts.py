"""Host names from PLSDB, PlasmidScope and GenBank/RefSeq, reduced to species and genus."""
from plasmidann import hosts


def test_a_binomial_gives_species_and_genus():
    assert hosts.normalise("Escherichia coli") == ("Escherichia coli", "Escherichia")
    # PLSDB writes underscores; a subspecies or strain is the same species.
    assert hosts.normalise("Klebsiella_pneumoniae") == ("Klebsiella pneumoniae", "Klebsiella")
    assert hosts.normalise("Salmonella enterica subsp. enterica serovar Typhimurium") == \
        ("Salmonella enterica", "Salmonella")


def test_a_genus_level_name_counts_only_toward_the_genus():
    assert hosts.normalise("Acidovorax sp. DW039") == ("", "Acidovorax")
    assert hosts.normalise("Candidatus Azobacteroides pseudotrichonymphae") == \
        ("Candidatus Azobacteroides pseudotrichonymphae", "Candidatus Azobacteroides")


def test_a_name_that_is_not_an_organism_is_no_host():
    for name in ("human gut metagenome", "uncultured bacterium", "Plasmid pDS56", "-", "",
                 "uncultured prokaryote", "bacterium"):
        assert hosts.normalise(name) == ("", ""), name


def test_a_higher_taxon_or_a_free_text_word_is_no_genus():
    """Names found in the analysis set: each would otherwise give a family, order, class or
    a non-taxon as the genus and inflate genus_count."""
    for name in ("Enterobacteriaceae bacterium", "Rhodobacteraceae_bacterium_SC52",
                 "Hyphomicrobiales bacterium 7MK25", "Mollicutes bacterium LVI A0006",
                 "Alphaproteobacteria bacterium AO1-B", "Nostocales cyanobacterium HT-58-2",
                 "Natrialbaceae archaeon AArc-T1-2", "Abditibacteriota_bacterium",
                 "Enterobacteriaceae endosymbiont of Donacia simplex",
                 "arsenite-oxidising bacterium NT-25", "endosymbiont_of_Sipalinus_gigas",
                 "glnQ allelic exchange vector pAG101",
                 "Sym plasmid", "Birmingham IncP-alpha plasmid", "Tomato aspermy virus",
                 "Affertcholeramvirus CTXphi", "Plautia stali symbiont",
                 "[Limnothrix rosea] IAM M-220"):
        assert hosts.normalise(name) == ("", ""), name
    # A genus before "endosymbiont" is the genus of an unnamed species, and a bracketed
    # misplaced genus is a genus.
    assert hosts.normalise("Wolbachia endosymbiont of Drosophila") == ("", "Wolbachia")
    assert hosts.normalise("Rickettsia_endosymbiont_of_Culicoides_impunctatus") == (
        "", "Rickettsia")
    assert hosts.normalise("[Clostridium] innocuum") == (
        "[Clostridium] innocuum", "[Clostridium]")


def test_a_phytoplasma_common_name_is_the_genus_candidatus_phytoplasma():
    """Names found in the analysis set: the first word is a plant, not a genus."""
    for name in ("Onion yellows phytoplasma", "Periwinkle little leaf phytoplasma",
                 "'Catharanthus_roseus'_aster_yellows_phytoplasma"):
        assert hosts.normalise(name) == ("", "Candidatus Phytoplasma"), name
    assert hosts.normalise("Candidatus Phytoplasma mali") == (
        "Candidatus Phytoplasma mali", "Candidatus Phytoplasma")


def test_a_list_of_hosts_gives_the_first_that_names_an_organism():
    """PlasmidScope joins the hosts of several source databases with commas."""
    assert hosts.normalise("Escherichia coli,Shigella flexneri") == (
        "Escherichia coli", "Escherichia")
    assert hosts.normalise("Klebsiella pneumoniae,-") == (
        "Klebsiella pneumoniae", "Klebsiella")
    assert hosts.normalise("-,Lactococcus lactis") == ("Lactococcus lactis", "Lactococcus")
    assert hosts.normalise("Salmonella enterica subsp. enterica serovar 4,[5],12:i:-") == (
        "Salmonella enterica", "Salmonella")


def test_the_first_source_that_names_an_organism_wins():
    assert hosts.resolve([("plsdb", ""), ("plasmidscope", "gut metagenome"),
                          ("organism", "Escherichia coli")]) == \
        ("Escherichia coli", "Escherichia", "organism")
    assert hosts.resolve([("plsdb", "Shigella_flexneri"),
                          ("organism", "Escherichia coli")])[2] == "plsdb"
    assert hosts.resolve([("plsdb", "")]) == ("", "", "")


def test_a_name_is_eukaryotic_only_when_every_taxid_it_names_is_in_eukaryota():
    """Against a small taxdump of the NCBI layout. "Bacillus" names a bacterial and an
    insect genus; "Bosea" plainly names a plant genus, and as "Bosea Das et al. 1996" a
    bacterial one; a name absent from the taxonomy is not eukaryotic."""
    import pathlib
    taxdump = pathlib.Path(__file__).parent / "data" / "taxdump"
    names = {"Saccharomyces cerevisiae", "Homo sapiens", "Escherichia coli", "Bacillus",
             "Bosea", "Nonexistens"}
    assert hosts.eukaryotic(names, taxdump) == {"Saccharomyces cerevisiae", "Homo sapiens"}
