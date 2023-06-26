# CACO

CACO: A core-attachment method with cross-species functional ortholog information to detect human protein complexes

## Datesets

- STRING: data/string_dataset.txt
- GO (from mouse): data/human_mus_go_reviewed_count.txt

## Reference

- CORUM: data/human_main_all.txt

## Usage of CACO based on GO terms
You can directly test your own dataset on our web-server: https://csuligroup.com/CACO
### Run demo:
```
python CACO.py
```

### Apply CACO on a specific dataset
```
python CACO.py --input ./data/string_dataset.txt --GOfile ./data/human_mus_go_reviewed_count.txt --Core_threshold 0.4 --output ./results/BioGRID_network.txt
```

### Parameters

> --input: PPI network
>
> --GOfile: GO file
>
> --Core_threshold: the threshold in generating Cores
>
> --output: the result file


## We also provide CACO when the GO terms are missing
### Apply CACO on a single PPI network
```
python CACO_top.py --input ./data/string_dataset.txt --Core_threshold 0.4 --output ./results/BioGRID_network.txt
```
### Parameters

> --input: PPI network
>
> --Core_threshold: the threshold in generating Cores
>
> --output: the result file

## Concat
Please feel free to contact us for any further questions.

- Wenkang Wang wangwk@csu.edu.cn
- Min Li limin@mail.csu.edu.cn

## References
CACO: A core-attachment method with cross-species functional ortholog information to detect human protein complexes.