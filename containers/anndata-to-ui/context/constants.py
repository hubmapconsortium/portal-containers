class AnnDataLayer(str, Enum):
    SPLICED = "spliced"
    UNSPLICED = "unspliced"
    SPLICED_UNSPLICED_SUM = "spliced_unspliced_sum"


class Assay(Enum):
    def __new__(
        cls,
        key: str,
        salmon_option: str,
        secondary_analysis_layer: AnnDataLayer,
        barcode_adj_performed: bool,
        barcode_adj_r1_r2: bool,
        keep_all_barcodes: bool,
    ):
        obj = object.__new__(cls)
        obj._value_ = key
        obj.salmon_option = salmon_option
        obj.secondary_analysis_layer = secondary_analysis_layer
        obj.barcode_adj_performed = barcode_adj_performed
        obj.barcode_adj_r1_r2 = barcode_adj_r1_r2
        obj.keep_all_barcodes = keep_all_barcodes
        return obj

    def __str__(self):
        return self.value

    CHROMIUM_V2 = "10x_v2", "--chromium", AnnDataLayer.SPLICED, False, False, False
    CHROMIUM_V3 = "10x", "--chromiumV3", AnnDataLayer.SPLICED, False, False, False
    SNARESEQ = "snareseq", "--snareseq", AnnDataLayer.SPLICED_UNSPLICED_SUM, True, False, True
    SCISEQ = "sciseq", "--sciseq", AnnDataLayer.SPLICED_UNSPLICED_SUM, True, True, True
    SLIDESEQ = "slideseq", "--slideseq", AnnDataLayer.SPLICED_UNSPLICED_SUM, True, False, False