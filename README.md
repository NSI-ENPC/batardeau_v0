# Code éléments-finis - Batardeau

## Présentation

Ce projet correspond à l'implémentation d'un code éléments-finis pour le calcul des batardeaux, à l'aide d'éléments de coques minces et de poutres (liernes, butons, etc etc). La structure est considérée encastrée ou en appui rotulé en pied.

Le moteur est codé en <code>C++</code>, et utilise la librairie **Eigen** sous licence MPL2. L'interface est codée en <code>Python</code> à l'aide des modules **Dash** et **Plotly**.

## Configuration VS Code

La librairie **Eigen** doit être présente sur le système, et son chemin connu. Le *compiler* utilisé sous Windows et Linux est <code>gcc</code>. Il convient de modifier les fichiers <code>tasks.json</code> et <code>c_cpp_properties.json</code> pour indiquer au *compiler* et à *Intellisense* où se situe la librairie mentionnée.
