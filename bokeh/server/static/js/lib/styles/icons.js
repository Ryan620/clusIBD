import "./root";
import * as _a from "../core/dom";
_a.styles.append(".bk-root .bk-tool-icon-box-select {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg0kduFrowAAAIdJREFUWMPtVtEKwCAI9KL//4e9DPZ3+wP3KgOjNZouFYI4C8q7s7DtB1lGIeMoRMRinCLXg/ML3EcFqpjjloOyZxRntxpwQ8HsgHYARKFAtSFrCg3TCdMFCE1BuuALEXJLjC4qENsFVXCESZw38/kWLOkC/K4PcOc/Hj03WkoDT3EaWW9egQul6CUbq90JTwAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-box-zoom {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg82t254aQAAAkBJREFUWMPN11+E1FEUB/DPTFn2qaeIpcSwr5NlUyJiKWVXWUqvlUh/iE3RY9mUekkPPURtLKNRrFJEeuphGfUUaVliiX1aVjGs6aG7+XX9ZnZ+d2fTl2vmnHvPPfeee/79Sk+may2/UQq/q7Qu+bAJoxjHIKqB/wlfUMcMVqI9bLZ+DGIKwzlzQ2GcxCx2xwvKOUKlaHTiX8bHNspjDONHkOmJBW5jIof/FvPh/06MZOb6cRc7cGn1AKUE5cdzlM/gAr5F/O24H3xkFRfxAbVygvK+cIsspjGWo1zgjeFpxL+BvnLw7laBA4xjIFJwrgu52DoVjKdY4HBEX8dSF3JLYe1fe6UcYCii3xWQjdfuSTnAtoheKCC7GNED5Zx4L4qt61jbTLHA94geKSC7P7ZeShQ0Inoi1IJuEOeORooFXkV0FZNdZs5qvFfKAeqYy7nZ6yg//HG0MBfffh71lFrQDCW2EvEP4mt4okZUDftz9rmGZkotmMxJRtlisy+MTniAWrty3AlXw0hFM2TD89l+oNsoOJXjbIs4EpqNtTCLXbiZ0g+M4mFObj8U3vsNjoZCVcmk60ZwthpepLZkB/AsivWfOJZxtpUQHfWib7KWDwzjeegBZJSdKFiE2qJTFFTwElsi/unQ/awXrU4WGMD7nOJxBY/1EO2iYConq93CHT1GOwucjdqnRyFz+VcHmMNefMY9nNkA3SWUOoXhQviSWQ4huLIRFlirFixnQq/XaKXUgg2xQNGv4V7x/RcW+AXPB3h7H1PaiQAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-zoom-in {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEgsUBmL8iQAAA2JJREFUWMO9l12IlFUYx3//MzPrLpSjkm5oN4FFIWVEl66IQlFYwtLOzozsjHdGRSCRF0sfBEXRVV0FQuQiLm5CZNBFgRRaRLVFhbJ2EdiN5gbK7toObTPn6eYdPTvNzPvOBz5Xh/ec5/n/n89zXtEHmZqeSXSuXBz/3zfdKvBWJHQrwZuRcP0El+QkbQXeBX6WZEgm6TtJk5lM5o4Lc+cV6qpf4Ga20Tm338zeATItVK9Ker6yvPzp4NDQ3+XieGsCU9MzTYumGbhz7m4ze9/MHgvBgItACrgfGAj2jgAvAYs3wlEujjc13kii8YyZrXXOfWhmo9GnFUlvOOemarVapVqtkslksmb2KjARqL62ecuWN9NxbRInzrldAXhV0uFSIfdew7G/gNLU9MwS8CwSmE3Oz88fcXG5blfpqVRq0Ix8VIAAX0XgrVL7HDCHGcCaWrV60LUBN8Dae58aQIxEqcA592I9M610JL0cpG/U9TIHJNKY3RV5z0R+7Nd4HZ0P1g/2RMBuegLAsRMnb4vT8d5vqKfMzOgtAlADrkmqGywmiMBTwfr3dC9j1Xv/r6Tvg/5/5ejxE6cO7M9faVbQZrYNOFSPmqQvVo9FKexvi5uWX58943aM7DwAfBDY+FbSCxP5sdkGx55GeguzrUEXPaSo2pFkAbiSZQCAzZJOmdkjwd6SpB/M7KykQTPbA2wDhoIzRzcNDx9MJwGNIXdJ0mEzmwbujL7dbma7gd03A7lKfnTOvf74nl0r6bonTUbujRSUCrm2d4L3/kvn3JPe+8+BDW2i9o+kT7z3kxP5sYsA6W47oE64TsR7P9tQL4vA2mh9WdIscKxUyJ0M7aR7acOGzikD65EQLEjaa2ZXzMwDFeB6qZBbbLTRE4EGeSaozNOZgYFf8qP7lmIvs354n0qlHpB0T7B9Ogl4IgJJrmjv/SiQjbrkD+BMUkfSbYATPdckrTOzkciWAXOlQu5cYgLdPEIapud9wMOR9zVJH3ViKx333mtHMJvNuoWFhZ3A+ojMcja77njXBEKwJJfTcqUyCIQ34Mf7nnh0paMnXacFuGoC1mr3AtuDfLzd8Zuyl+rfuGn4HLAD+Az4qZQf+61TAj0Noj8vX6oC35SL43u7teG6rf5+iXppwW7/JUL5D03qaFRvvUe+AAAAAElFTkSuQmCC\");\n}\n.bk-root .bk-tool-icon-zoom-out {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEgsHgty9VwAAA0FJREFUWMO9l09oXFUUxn/fmXlpItppi22k7UJBRSlVkCytSAuKUloIdjKT0El3FXVXdVFKRVAQV7qQohsNwdA0UFvBhYtqUVyIVlRaogtFQVq7qSTVjA3z3nHzBq/jvPmTN/Ss7rv3nvN99/y794kByMzcfE/7picn/jenmwWeRUI3E7wdCRskuCSTdDfwBvCtJEdySV9KOhpF0e0/LF5SqKtBgbv7ZjObcvfXgShD9Zqk5+orKx8Oj4z8NT05kU1gZm6+bdK0Azezu9z9hLs/HoIBvwAF4H5gKFh7B3gBWFY3460kWve4+3oze9fdx9OpVUmvmNlMHMf1RqNBFEUldz8OHAxUX9q6bduryut+Sfvc/Wz62ZD0fK1afjND9y3gGSRwv1GMojstTxUUCoVhdyopEYDzKXjWwZ4FFnEHWBc3Goet00m7lZlZYQixKw0FZnakGZksHUnHgvCN5/KARBH37enpOVg58H13HV0Kxg/kIuD/ngSA2ZMLt3bTSZJkUzNk7k4+D0AM/CGpaXCyBw/sC8Y/qZd2GpZiuL9YLN4Sx/HpoP5/c/exQ1OVq+1yyt13SLoArEsJnMjlgfOffvK3u58Kprab2QezJxfG2iTzUzI70wRPG9jbmpmb95SNB9mpzp7/j2yVdNbdx4K565K+cvfPJQ27+x5gBzAS7Hlvy+jo4WIvoC3kWpcvS3rR3eeAO9K529x9N7C7zX6AC2b28hN7Hl1Vt44niVq13LUjmtlYkiQfA5s6eO+GpDNJkhw9NFX5ueNt2ARodyF1IHIN2JiOl4H16fiKpK+B2Vq1vBAqFAf4IJkGNiIhWJK0192vunsC1IE/a9XycquNXARa5OnApeeioaHvKuP7r3dTGsiLqFAo7JR0T7B8rhfwXARa2us4UEqr5Ffgs151i/08oTNKdIO770ptObBYq5Yv5ibQq/sl3Qc8lJ4+lnSqH1vFfp9koZRKJVtaWnqkWXqSVkqlDe+vmUDWpZMlK/X6MBDegKf3P/nYaj8ErN9fqZBYEsf3Ag8G8Xit33BaniTcvGX0IvAw8BHwTa1y4Md+CeRqRL9fudwAvpienNi7Vhu21uwflOT+L+i1X2TJP57iUvUFtHWsAAAAAElFTkSuQmCC\");\n}\n.bk-root .bk-tool-icon-help {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAAXNSR0IArs4c6QAAAAlwSFlzAAALEwAACxMBAJqcGAAABltpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6dGlmZj0iaHR0cDovL25zLmFkb2JlLmNvbS90aWZmLzEuMC8iCiAgICAgICAgICAgIHhtbG5zOmV4aWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20vZXhpZi8xLjAvIgogICAgICAgICAgICB4bWxuczp4bXBNTT0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL21tLyIKICAgICAgICAgICAgeG1sbnM6c3RSZWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9zVHlwZS9SZXNvdXJjZVJlZiMiCiAgICAgICAgICAgIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIKICAgICAgICAgICAgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIj4KICAgICAgICAgPHRpZmY6UmVzb2x1dGlvblVuaXQ+MjwvdGlmZjpSZXNvbHV0aW9uVW5pdD4KICAgICAgICAgPHRpZmY6Q29tcHJlc3Npb24+NTwvdGlmZjpDb21wcmVzc2lvbj4KICAgICAgICAgPHRpZmY6WFJlc29sdXRpb24+NzI8L3RpZmY6WFJlc29sdXRpb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOllSZXNvbHV0aW9uPjcyPC90aWZmOllSZXNvbHV0aW9uPgogICAgICAgICA8ZXhpZjpQaXhlbFlEaW1lbnNpb24+MzI8L2V4aWY6UGl4ZWxZRGltZW5zaW9uPgogICAgICAgICA8ZXhpZjpQaXhlbFhEaW1lbnNpb24+MzI8L2V4aWY6UGl4ZWxYRGltZW5zaW9uPgogICAgICAgICA8ZXhpZjpDb2xvclNwYWNlPjE8L2V4aWY6Q29sb3JTcGFjZT4KICAgICAgICAgPHhtcE1NOkluc3RhbmNlSUQ+eG1wLmlpZDpBODVDNDBDMzIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRTwveG1wTU06SW5zdGFuY2VJRD4KICAgICAgICAgPHhtcE1NOkRvY3VtZW50SUQ+eG1wLmRpZDpBODVDNDBDNDIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRTwveG1wTU06RG9jdW1lbnRJRD4KICAgICAgICAgPHhtcE1NOkRlcml2ZWRGcm9tIHJkZjpwYXJzZVR5cGU9IlJlc291cmNlIj4KICAgICAgICAgICAgPHN0UmVmOmluc3RhbmNlSUQ+eG1wLmlpZDpBODVDNDBDMTIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRTwvc3RSZWY6aW5zdGFuY2VJRD4KICAgICAgICAgICAgPHN0UmVmOmRvY3VtZW50SUQ+eG1wLmRpZDpBODVDNDBDMjIwQjMxMUU0ODREQUYzNzM5QTM2MjBCRTwvc3RSZWY6ZG9jdW1lbnRJRD4KICAgICAgICAgPC94bXBNTTpEZXJpdmVkRnJvbT4KICAgICAgICAgPGRjOnN1YmplY3Q+CiAgICAgICAgICAgIDxyZGY6U2VxLz4KICAgICAgICAgPC9kYzpzdWJqZWN0PgogICAgICAgICA8eG1wOk1vZGlmeURhdGU+MjAxNjoxMToyOCAxMToxMTo4MjwveG1wOk1vZGlmeURhdGU+CiAgICAgICAgIDx4bXA6Q3JlYXRvclRvb2w+UGl4ZWxtYXRvciAzLjY8L3htcDpDcmVhdG9yVG9vbD4KICAgICAgPC9yZGY6RGVzY3JpcHRpb24+CiAgIDwvcmRmOlJERj4KPC94OnhtcG1ldGE+Cphjt2AAAAT7SURBVFgJxRdbaFxFdGb2bhui227BWrsVKYgf2kJUbP9EUPuzEB803WTXJjH61Q/7Ya1+CMYKEVTsh4J/EpvY7BoabUiNiA8s1p+4KIhpoUUEselHqyS76TbZ3HuP58ydc3d2u4+IkQxczpz3mZkzZ86VYpXjvenpjZsLhUcliE4AuUuASAgptmt1EFdwPiclzIIUUwubNn17OJlcXo1p2UpodHRiux9xB1Eug1+slbzhFxGOKc851tu7/0oznYYBDA8Pt0U2tL8KQryIq2tvZqQhD0QJHRz3yqWhgYGBpXpydQMwqz6NCnurleCSADkJEfgKfOePqL80R/wV1ZaQyr1LenKfkPCkEPKeaj0xg7vxVL3duCmA0Vyuw/fl52hgBxsBED+h4Cv9z3R/zbRm8MTJTx7HQN7GQB6w5C4L4SX7M5lfLBpurjXMyvNIShiyi0l1pL8n9b7EDGPR8fHxzSsQ6XDB3618/xqo6Pk25V5MpVJllgHM1BO58RdQ612kOYZ+GXdij70TYQB05mpj+1kU5G2fB+l3PZtOf8NGx6ambnMXb3yAxg8wjSEG6OKKR9oicBQD+ZvpH2Wzj0lQpxCPG9qMv1x6hHNCsSAlHM7ZOa682vlI9tRDbvHGbD3nZAPpDoD/3JIrLpAs26UFkC3EMUA99hpfGtEBfJjNJnS2Gwnadnvl+Xw+iuc3DAJuNyIaSCHpilVldyDjjUxj3WDZIAhxhHHyRcdNuA7AAfUaXzVKODpzFiZ4/uLvh5G+m2no+C/pyIf7MqlEJB7bpqR6nXkEUfbeawuLaZsW2ISfNQ2vtaktQlGFQyIVGT0o2+2EC4iQNGwjBIN9qdQ5Qg4mk4X4rW3vCClLtowE2FOFUxKDfNmiZci3ovKKRFPh4FK9q4Zbdr+lKKJiA13TcHR2dmLBgdmQ0GAS2MZaEowY+XbAk09IvgtYZGp16SyvFhaHcIUh645t8T9DBCcnz5zZ4hZLu3DzK2QlL1QQa0Y+pHiJKPSuOGj3PmZTheM5w2TwqBxnvBZOTk7G5gvXJ5Aelms8wnJURL+olSWcfEhf6gDoUXPMq6ZlqbzWU2pE+3hi4s6F68tfIj9cBMlikr7Z0/P0b/X0yIcUXsDCF1WhtL4OROHaXk+xlkbV0Cu732Nmhc4peaWSg73pA8dq5RkvO37ldUTfXCKZv2q45MkhvG87WQEzpCCUSvV1d9GONBy3lMvgKSwrZig8gjAietWY0QriylO2jIo4yVbOSb7KB/qmI9BPKjHpSSXYauRyn92Nq9/Kcrj13x3s3v8D481glQ/0raiNYgX9njPSBOImbrHZePl+tfFmc9sH+Xaoh8NjOKSVdDMhjjYzQLy+dFceH5+IJQf9VYXX4tROg4ZFU8m31M3mfPEqUoJqCGJfvWpo2xnNfdrhC28n06SCeSzNZxlvBINGRXCtKS7EY1uV6V7HWAm38y1cXaXsMcOCvr9ySPj+af7A1U2HJXHzVNvUXVLIGyPf+jV0pf8GHoN+TLAyPkidTCi2RpPApmnR0Bd1zGRaB/B8Oj2HSw7LLbVR1MmskW8RdEWVXSJf3JbpAMgRtc4IZoxTh9qotQjCasm46M0YX9pV1VmbpvRH5OwwgdRtSg2vKaAz/1dNKVtb17Y8DCL4HVufHxMOYl1/zTgIgiYvBnFKfaNp3YjTdPz3n9Na8//X7/k/O1tdwopcZlcAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-hover {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4oVHp0SwAAAQJJREFUWMPtlsENgzAMRb8RQ5VJItFDOgaZAMaAA0iZpN3KPZSoEEHSQBCViI/G8pfNt/KAFFcPshPdoAGgZkYVVYjQAFCyFLN8tlAbXRwAxp61nc9XCkGERpZCxRDvBl0zoxp7K98GAACxxH29srNNmPsK2l7zHoHHXZDr+/9vwDfB3kgeSB5IHkgeOH0DmesJjSXi6pUvkYt5u9teVy6aWREDM0D0BRvmGRV5N6DsQkMzI64FidtI5t3AOKWaFhuioY8dlYf9TO1PREUh/9HVeAqzIThHgWZ6MuNmC1jiL1mK4pAzlKUojEmNsxcmL0J60tazWjLZFpClPbd9BMJfL95145YajN5RHQAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-crosshair {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAADEUlEQVRYR81XXVIaQRCeHqug8CXmBNETaE4gniDwIgpVspxAbxC9ATkBkCpQ8gKeQDiB5AQxNyAvUlrldr7eHxyGXZi1rMJ5opbp7m++7un+htSGF204vsoMoNXrlzSpfWa1oxQfhAegCZGaEtPorHo8znIoJwCt6+td8uk7ApUQCIHTF4BNAWzImq8ap6cP68CsBdDp9i9ZqXM7ML79g/EnCWD+jgMKENKqWT+tXK0CkQqgNRjs0OxpQIqKhoMxaG6/6JeRnK7T6yO2UvVqhYSlLX+ryORfgKn9ORDFIy7ky41yGcwsr0QAQfDH5zucOswx819fs4egI9OFCcD8DjBF7VNbEX0JzdWEt3NHSSASAcCxBDqMgt/623kvyTgNgNjJIfTjk4D4FqaJR1715MjmYAmA5Bx3AwUXQL+t105KaTlcBSC26XRvhjEIoLiq1yqXpr8FAGG16/ug4IT27fxBWu7EiQuAiImJpEMKE6nYM30uAIDDttSUOPfJP7JzbjPhAiBIh9QE67vIvoOi9WJfCwDavf40ulpjbCqmUf+W753ezURuh7Dg1SqflwAEHU6pgfyBq9Y4qx0LG++2fnZ/eUzcstmdM2AWH+jfc+liWdBJfSENf8Lifi3GVwC9mybOfi5dzatWVrbbLIHNva8p5h/16gkaFiLGGxbufkoE6XguwePiXLF3XmMfCUCUAqtKXU7sumd1CowOuJEi3Pg1FBpjitIGhyvVSfvmjci6ZR+rFQfDiPVE2jFYeICQ+PoewwjC5h7CZld6DBdyu6nDSKgzOyIMhmhK5TTqXYbRorZYM46TmpKAAOrGWwSJJekSB1yqJNOzp1Gs7YJ0EDeySDIMtJbQHh6Kf/uFfNFZkolJICRmz0P8DKWZuIG2g1hpok+Mk0Qphs0h9lzMtWRoNvYLuVImUWrmPJDlBKeRBDfATGOpHkhw670QSHWGLLckmF1PTsMlYqMJpyUbiO0weiMMceqLVTcotnMCYAYJJbcuQrVgZFP0NOOJYpr62pf3AmrHfWUG4O7abefGAfwH7EXSMJafOlYAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-lasso-select {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEgwlGP1qdAAABMBJREFUWMO9V1uIVVUY/r61z57ZMx4DnbzgkbQXL5iCJphlWdpIGY4jpFBkEiU9ZNaDRRcITcIwMwgxoQtU2IMXdAZfMjFvpERXYiSbysyBEXFmyuHMnLP32uvrwT2xnY5nxvHQ93Jg7fWv71/r//7L4a59TRgqJk+Z6v3a+sv0OI5nk5wu6VaSVZImAThHsgjgrKTvM5nMUWvtmf5n8HodCIKgOgzDhc65pSTrJQWDsSNpJX1ljHnDOfdT37oZLLHv+8OMMasKhcIJ59xHAJYMlhwAJGUAzJfUTHLFuFzOG5QDU6dNMyQfs9Yedc5tBpAD4IYYNQGoBrDtQnt7/b0LFrJsCHzfn2itfQfAnZLiazytA3AaQAuAiwDaEgeNpGkkswAWSBqRONB38b88z5uTKePt6iiKXkk8jq+iJC5LOmiMaTLGHLPWhmWeHr7vV0dRtATAapAzIVmSo51zyzIlbm2stesFPA6pKk0r6Ryg93y/ek8YFvPOOTg3cDSiKCoC2OP7/rEoirYm4rUkF12lAWNM1lr7lqQn0+QA8gI2jBg5cj6Aj8OwmB+KAKIoukhyp6SRJAUgl0ndPLDWPi9pJQCbuviXvu+/GIZhW1dnJ24UJFuTjCCA2ADA8sYGWmsXS3qmL94kDYAtkh4Nw7ANlQJ5U6INT1KrAYC9zQdykl7nFSj5fXp5Y8NWVBhy7mUAjqShMYdMXV2dJ2klyRwAJ8lIeuGWCRMP7N7frEqSG2OmAFhKshNAp5wrmO7u7jEAngPQm1S2z2pqapr+OPt7XEly0oxwzq2RdFmSD2AMgKKJouhhAL4kA+Cs53l7e3t7uytJHgRBreTWkXwkKVJnJD0B4GAGwIJE9R6AFufc6UqSZ7PZbD6ff5dkA4CQZEHSqwAOISmXtwGIE+F1SeqqIP8d+Xz+C0mLJYWSAODteXffczjdDQNJ0BWMCoLg5gqIbRTJNwHsljQhUb0luWPM2LE7Thw/9m/5NCT/TByxAOYWi8X6/gdWV1dnfN8fNRBxJpMZTXKdc+6IpFVJWAEgkvSJpA0X2tvtVTaSjgOYBCAEEADYSHK87/sfhmEYA9gShuEDkgzJHyWtB/B1irQ2juP7ADxkrX0wOUOpzmdpzEY590HJ7Ni1r2kSyZOSiv2+hSRjSTXp/QAukzySNJOJkmalyNIl10hqMcasdc61XDNcQRD8BnITgNp+36r6kfcNFMMlLQGwTNLMEuQGQBfJl2bdPru+HDkAZAqFQux53jZHEsC6aw0eg2gylNRBcqcx5v04ji999+03AwsWAOI4Lsy9a94WkisAnE5a5WCJYwCfA1g7LJudI2lTHMeXBm1faiQzxkyRtF3S5CTupeAB+KG2tnZFT0/P30NO2VKLzrmfAbwGMipjG5Oc0dPTc0Md05SZ5U4Q2FxChErtEYD7jTGNQ3UgM8Asv90Yc9I5LSKRlXSI5CxJa0jWSALJjKRnAewfkniT+vwf7N7fXHK9rq7O7+jo+BTA/NRrdBpjnnLOnUrvXd7YMPQXSBunneno6IhIHgYwW1JtkgmBpBkATlVMAwOk3nFJ+VSoqgCMr6gIy2FcLtdKspAedyQN/98caDt/3kpyabUmf8WvG/8A1vODTBVE/0MAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-pan {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4lKssI9gAAAOtJREFUWMPVll0KwyAMgNPgoc0JzDX2Mtgp3csKErSamGabIEUo/T6bHz0ezxdsjPJ5kvUDaROem7VJAp3gufkbtwtI+JYEOsHNEugIN0mgM1wtsVoF1MnyKtZHZBW4DVxoMh6jaAW0MTfnBAbALyUwCD6UwEB4VyJN4FXx4aqUAACgFLjzrsRP9AECAP4Cm88QtJeJrGivdeNdPpko+j1H7XzUB+6WYHmo4eDk4wj41XFMEfBZGXpK0F/eB+QhVcXslVo7i6eANjF5NYSojCN7wi05MJNgbfKiMaPZA75TBVKCrWWbnGrb3DPePZ9Bcbe/QecAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-xpan {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4X4hxZdgAAAMpJREFUWMPtlsEKwjAMhr/pwOOedINJe/PobWXCfAIvgo/nA4heOiilZQqN2yE5lpD/I38SWt3uD9aMHSuHAiiAAmwaYCqoM/0KMABtQYDW11wEaHyiEei28bWb8LGOkk5C4iEEgE11YBQWDyHGuAMD0CeS30IQPfACbC3o+Vd2bOIOWMCtoO1mC+ap3CfmoCokFs/SZd6E0ILjnzrhvFbyEJ2FIZzXyB6iZ3AkjITn8WOdSbbAoaD4NSW+tIZdQYBOPyQKoAAKkIsPv0se4A/1UC0AAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-ypan {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4anK0lywAAAMVJREFUWMPtlzEKwzAMRX/S7rlpIMXeOnaLaME36FLo8XqCdNFghGljyc4kgQi2Q/SUj0F/eL7eMMTKz6j9wNlYPGRrFcSoLH4XxQPvdQeYuPOlcLbw2dRTgqvoXEaolWM0aP4LYm0NkHYWzyFSSwlmzjw2sR6OvAXNwgEcwAEcwAEcwAEcoGYk20SiMCHlmVoCzACoojEqjHBmCeJOCOo1lgPA7Q8E8TvdjMmHuzsV3NFD4w+1t+Ai/gTx3qHuOFqdMQB8ASMwJX0IEHOeAAAAAElFTkSuQmCC\");\n}\n.bk-root .bk-tool-icon-range {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAAXNSR0IArs4c6QAAAAlwSFlzAAALEwAACxMBAJqcGAAABCJpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6dGlmZj0iaHR0cDovL25zLmFkb2JlLmNvbS90aWZmLzEuMC8iCiAgICAgICAgICAgIHhtbG5zOmV4aWY9Imh0dHA6Ly9ucy5hZG9iZS5jb20vZXhpZi8xLjAvIgogICAgICAgICAgICB4bWxuczpkYz0iaHR0cDovL3B1cmwub3JnL2RjL2VsZW1lbnRzLzEuMS8iCiAgICAgICAgICAgIHhtbG5zOnhtcD0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjU8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOlhSZXNvbHV0aW9uPjcyPC90aWZmOlhSZXNvbHV0aW9uPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICAgICA8dGlmZjpZUmVzb2x1dGlvbj43MjwvdGlmZjpZUmVzb2x1dGlvbj4KICAgICAgICAgPGV4aWY6UGl4ZWxYRGltZW5zaW9uPjMyPC9leGlmOlBpeGVsWERpbWVuc2lvbj4KICAgICAgICAgPGV4aWY6Q29sb3JTcGFjZT4xPC9leGlmOkNvbG9yU3BhY2U+CiAgICAgICAgIDxleGlmOlBpeGVsWURpbWVuc2lvbj4zMjwvZXhpZjpQaXhlbFlEaW1lbnNpb24+CiAgICAgICAgIDxkYzpzdWJqZWN0PgogICAgICAgICAgICA8cmRmOkJhZy8+CiAgICAgICAgIDwvZGM6c3ViamVjdD4KICAgICAgICAgPHhtcDpNb2RpZnlEYXRlPjIwMTgtMDQtMjhUMTQ6MDQ6NDk8L3htcDpNb2RpZnlEYXRlPgogICAgICAgICA8eG1wOkNyZWF0b3JUb29sPlBpeGVsbWF0b3IgMy43PC94bXA6Q3JlYXRvclRvb2w+CiAgICAgIDwvcmRmOkRlc2NyaXB0aW9uPgogICA8L3JkZjpSREY+CjwveDp4bXBtZXRhPgrsrWBhAAAD60lEQVRYCcVWv2scRxSemZ097SHbSeWkcYwwclDhzr1Q5T6QE1LghP6BGNIYJGRWNlaZItiFK1mr+JAu4HQu0kjpU8sgF3ITAsaFg0hOvt2Zyfvmdsa7a610Unx44Zgf773vvfneezPHNzrbhn3CT3xC3wPXYOC8LDzqdi8YY/gwh4BeknS/2th6dr2kf94AOp3OFyWgMyziOPbMDxV9FTtJnl1ut795Xd0/YQ0/vtYQwMT1KXWCfr2IjOWwtNehwN4xL9ykTrm6Pzl58yLn3J+mKh9mXbT3uRjGEDph+O8/TjfP5dBp7Ha7AX7O3o5nZeD/0E/OGyXntDgzA0X6qmCnrVutVlrUWV9f/3xo+pwhGDhvEPHOjoxnZjJggXmMHzBQ7NGNp9vxk61fr0HR7e/u7pZzCGHlc7qwBYYTT7tJYSx1AQzppyFPft5apta9w7SKcn0b7P7+/jCsDQ5mbc0dCmIJGDN0ehdcjsmkm6A6KUeKFOTE11PLxrC7Ukqh3ylL2fT0NAP9q6ur6rRCJJYsbKB0JsbCKMuy+xREePDyxQPCz+Crlw062QcA5wBOOt1l6vIl2WiI9F1fN6Q+BBqit6hEC4Hk08GQJMn4myjSP7RavVxgdaVUh/3U6HCMsPr9pYnJKRziHtWQ+un58+hGs6nsjQSjpuTyKGN3CX+FBwHXSiEVgjP+O8X6N12kIePES+GzTKAkGbNp8yJsGUMVzz8jPKReiyAQRimy5/cjye5RpF8utFp/+nwmT7d/NMzcFkS7yjJNGDaPURQxIQThEQy0SyF4l5WJYYhBa816vZ6dU7A6CAhbZVow/pDe0O9hVOoCi13r4BgBAvJHqMSQL2vE/iH6IAXEwgrRVUmBoRRwnwJQT98xEeVeSUyB4dJ5nwJBKdCFFGRmUCcu7rwIYypCTblaChuNBhWODrman5ub+4v0rMNBt8z6Ezh7GksJQpCbm79cMQE7QBFm/X6f0rjWnv8WRYg/QdbUpwDAEBy8vPyA8rNGzg3a8MiElwiM7dAtRqNoNptjGPM1laVxP9umWEMGLOKhKUOJDtBwDmzsw9fC/CzHr9SGuCTi2LbbKvVtmqXpCjMihBFa79Wrt5fGx9PDzc3fmu32Lf8qFliwU9emKhBSp+kRKn/hu9k1COEDbFdt/BoKWOAkuEbdVYyoIXv8+I/QK9dMHEb1Knb7MHOv8LFFOsjzCVHWOD7Ltn+MXCRF4729vWMDK+p8rLkvwjLg4N4v741m5YuwCI9CvHp1Ha8gFdBoPnQAkGsYYGxxcfEI7QQlFCTGUXwjAz4tWF+EpymOWu7fglE7qsOvrYE6g4+9/x/vhRbMdLOCFgAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-polygon-select {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEjc1OfiVKAAAAe1JREFUWMPt1r9rU1EUB/DPK0XbqphFHETo4OCiFhwF0V1KHbRSROLqon+AUMVRRFBwEbRFMBiV+mMW/wIxi5OD1kERRVKRJHUwLvfBTZrU5OWBGXLgQu7Jfe98z/ec7z0vKa88b2q1BDtRHdAPBaylm1NzsxsOjPnPNt6WSWprbft+/c3I3zOAjhT1Y4+fvcjEQJIXnVECSa+AhqIHqlHH5lWCZoe+Gk4GRgDG86j9SAUdlDBSQaZhlOkuHyoVdJmsw98D1S5fM4NYM1LCpqM+Lwa240oLgmZzpVZvzKT75VLZcqksSZKWlQeAy/iORVwIvh31xvotvK7VG3Px4aWHj3Jl4C2uYSvq+Bn8v6LLbaVWb9zsBiKLCvbiNG7gLm7jAYqbPHMJMziZ9lsKoh8GtqCEVVzHftwJn+TFHp4/hg8BSCYVfMOZoPEv2NZGdy9WCGUr9toDR3E2/H4V6nwRe/BmgN65H1ZhvMuB3XiKIyFoGefwO6ysVkUlrNUNsyAK/jli533Q+Y8cJFvAeXyMS1CI/jiMr/gUtD2LQwMGr4R3p7bY3oQHQ5b38CT4D2AXXg6YcQXHpyYnlqKsi5iOAVSwL9zd7zJ09r+Cpwq72omFMazjT9Dnibym0dTkRDUKrrgwH7MwXVyYB38BstaGDfLUTsgAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-redo {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4itK+dVQAAAaFJREFUWMPt1L1rFFEUBfDfJDaBBSslIFjbaSFp1FJQFMVCHkzhKIqdUYOCoBgErVz8rCwiTDMwBCIKipDWyip/gxAIWAmBgBC0eYFh2Gx2l9lFcA5M8e59782Zc84dWrT435Hs1siLchqn43MS0zgW22vYxjesYjVLw3YjBPKinMUTBOwf8J5fKLGYpWFjJAJ5Uc7gIW6jM6Kim3iNZ1katgYmEL/6I+YasvY7Lg6iRpIX5VF8wuEe/XV8wGf8jN6LWTiAc7iEQ7ucPZ+lYW0vAtfwvlbfwCKW9gpXDOv1mJvZHiSO91MiyYsyiQSuxtpXXM7SsDmM5nlRdrCMMz3sOJWl4Xevc/vwBzdwAl+yNNwZxfRI+GxelK9ikHcwh8d4NNR/YFRES1ZwoTYdR7I0rNf3TzVNIGbmSvR/Bx08mIgCFSVu4l2ltIWD9WxNGR+W8KOynqnZ0rwCeVG+wa0hjrxtWoF5dAfc28V8Mib/n+Nev5dnabg/zgw87aNEN/bHOwVRiRe4Wym9zNKwMKkpgIWKEt24njxiJlq0aPFv4i9ZWXMSPPhE/QAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-reset {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4gWqH8eQAABLdJREFUWMPtlktsVGUUx3/nfvfOlLQaY2IiRRMQIRpI0PjamJhoVASDvNpCpYw1vJQYSVwZwIVQF6wwRHmkAUof9ElrI6VqDAXcID4TF0IiYQMkSlTokNCZ+b7jove2t+NMH7rQBWd3v+989/zP+Z8X3Jb/WGQySvUNTQBJESkNguAVYIWqzhaRhwBU9WcR+QXoymazn6jqzUQiMQSQzWZRVdal1vwzAI2tHQBPOuc2AbWTdOyQ53n7nHNfRwee51GzqoIQMCLDpr3x/tLQ0oZzrk5Vj0/BOEBt+KYuOlBVGlrahr0Wob27t3gEjnZ2AyQzmUwHsDgP6J/AYRE553neDwDOuUdU9QngNeCumK4TkRMhZUORcYC1qysLA6iuSQHIwkWLD6lqapQsuSmwTVV3h99I7EcAR462A2xR2Ilq6ehTaejvO1774kuLNALR33eclsaGsQDe3fYegHl43vyNwEeqGl1963mm2jl7YZRTQ82qlWP4HM6ZToC5ztkW4LHQoALru7s6Di5dvlIj/e6ujrEAWoZDn8hmMjXATMACGaAVuBjXTVVXFc/AxhaA+4zvn1DV+eHxVWPMAmvtb5GeMWZyZVhI2rt7qVy2pOh9U1snwIPW2vMi4oWJuBPYHkVAVScPoKmtkzVVK6cEMsyJraHhiCqJqJUwj/JRz7TW1iSSyR2rVyylqa0Ta+24Ic8vXaAEmDFc/l5Z2A/80OibuVyuz/f9ElUdHCmvw82t5HK5h6y1PYhsz2YyGw43t2KtBZHIGwB6+j4rCkBVUdV7gXrggnPuu8h4eP+xMeZS2D0rJYZ6AdAMzAt1b4nI26p6IFZOY8pugijcKSIHVLUK0LyST4vnrVfnWr3mjmP4QTATaERkXkypRFX3isjmuHdRJEK6Ckqquopp06bdKCkp2Sgi7XnGLcg7gzeutwNIiPYc8HixqIrIOlU9ONVIhHPEd851icgSVXUiskVV94gIqoonIt0i8gfQCfwae38e6BWRXuBZz5jZ8VbaOE4EIqlZVUEQBLlkMplS1QER2RwkEnsSyaREDUzyeNsvIhvCMqkH1kdIJ2o+k8iJB1LVVRfjZ6nqqlEAIbdVQGto8Lrv+/dbawcjAL7vc+6bs+zetetfLSHxniIFGofGGsU2oC7eOCbDfZ7nQawBOSAX74SF9oEPImOq+r7nmVmxb5raukZa8UReGmNmhbMkAwwBH467EYVZe49z7kdgenj8k7V2oTHm8kgdWcvrNdVFjR8cHkYzjDH9wLjDaEwEzpwa4MypgWvAjtjxfGNMj4jMiT+M+kFsZI/Q6Pv+HGNMT8w4wI7TAyevxXVPD5z8+zD64tRXAMHVK1eaVLUyVvuDqroV2BOnJF4ZIedviUidqt4Re9s+vbx8zZXLl7PR2+nl5Tz/zNOFp2FzxzGAklw22wUsLLaSKXwf8vhosZUM6PeDYEUum70VHfpBwKsVyyfeikOP6oBNwN1TrLbfgX3A1kKLzKeff8nLLzw38T5wZDgxn1LnNk5lLRfP26/OnR2hwfNYW2Atn9RCsrf+EECyrKysDFimqhXhyjY3VLkAXBKRDqA7nU6nS0tLhyIj6XSaN9bVclv+l/IXAmkwvZc+jNUAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-save {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4UexUIzAAAAIRJREFUWMNjXLhs5X+GAQRMDAMMWJDYjGhyf7CoIQf8x2H+f0KGM9M7BBio5FNcITo408CoA0YdQM1cwEhtB/ylgqMkCJmFLwrOQguj/xTg50hmkeyARAYGhlNUCIXjDAwM0eREwTUGBgbz0Ww46oBRB4w6YNQBow4YdcCIahP+H5EhAAAH2R8hH3Rg0QAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-tap-select {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAA2hpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADw/eHBhY2tldCBiZWdpbj0i77u/IiBpZD0iVzVNME1wQ2VoaUh6cmVTek5UY3prYzlkIj8+IDx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1NzowMSAgICAgICAgIj4gPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4gPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIgeG1sbnM6eG1wTU09Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC9tbS8iIHhtbG5zOnN0UmVmPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvc1R5cGUvUmVzb3VyY2VSZWYjIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtcE1NOk9yaWdpbmFsRG9jdW1lbnRJRD0ieG1wLmRpZDo3NzIwRUFGMDYyMjE2ODExOTdBNUNBNjVEQTY5OTRDRSIgeG1wTU06RG9jdW1lbnRJRD0ieG1wLmRpZDpCOTJBQzE0RDQ0RDUxMUU0QTE0ODk2NTE1M0M0MkZENCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDpCOTJBQzE0QzQ0RDUxMUU0QTE0ODk2NTE1M0M0MkZENCIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ1M1LjEgTWFjaW50b3NoIj4gPHhtcE1NOkRlcml2ZWRGcm9tIHN0UmVmOmluc3RhbmNlSUQ9InhtcC5paWQ6OTQ0QzIwMUM1RjIxNjgxMUE3QkFFMzhGRjc2NTI3MjgiIHN0UmVmOmRvY3VtZW50SUQ9InhtcC5kaWQ6NzcyMEVBRjA2MjIxNjgxMTk3QTVDQTY1REE2OTk0Q0UiLz4gPC9yZGY6RGVzY3JpcHRpb24+IDwvcmRmOlJERj4gPC94OnhtcG1ldGE+IDw/eHBhY2tldCBlbmQ9InIiPz6eYZ88AAADLklEQVR42rSXf2TUYRzHv7tuGcfE6Vwb5zLSSjEj7Y9KWqfEmFZJP+yPMdKKmUrrn0iUfjhWlLFi6YfNrF+StBoTo39iYkTGco4xxxG59P7k/T2PT8/37nu3bx9ezvPj+zyf5/PreS78bGLS8SmrwE6yje3NHJsDBTALpknBz6JhH3NiYAB0gHqPOVv52wJ6QQ48BzdAttTioRJjdeA8mAHHS2xuk3p+M8M16ipVQE49Ds6CiFO9RLjGONf05QLx6wPQaBlbBlPgJVgkP0ETiIJ2sB/E1XfimjfgBOOlKDUqCGOcqBcQnw6BYW5YTo4wbvQhMmCfGRemC2rBiGXzWUb+kM/NRZ6CHWBM9ce5R61NgX6ayhSJ5EPlItlDRNkz4JbFHf06BkSzHjXxM+gDv1S/mPUo2AXWgt9UUHL/IVhS8yUV1/EbV3o4N+NaoE9Fu/i827K5pNYHnqAVJECShWmAaddpscYFFXwR7vnXBRGlnUN/L6kqKJlxnRUuDbaDBiL+vst5d4gpcpBrqk/2jIgCKVUolhntplzivHmwh4stGOPfwBWwl/2dpp8p7xjQZqFLiQJtauKkivYm+kzccpK57yXfOUe+P23JqAnVbhMFmlXntCWnxbT31am9ZJ4BJifsUmNTqt0cYhA5ypympPg7VkEKunPbVb8cIG+0kyHLJZNR7fUMooUKFHAPkfQo58VLK+RzwRDd4FdWG9mjpaAXzqkJa1R7kQttqEABWXMjOOxxVRfnhRm5URX1prk/0pQHwNcKlchZ+jdpC+hFdVqO0my9Hj5dkYgCn1Rfh/KdlNDHrJhPqlDih+IfBd6qwpOgEqYMsorJ2HtWxtagLJDn/W3KRfPOZhoeBJfZPgVeGKeKrkQBh5dLXl25Ny3pc4/1fkTdbvFqFQgbxWeYD0hXulhQ0pYiM1jG547fcbMQpVnHTZEn9W3ljsCzwHxCdVteNHIZvQa7/7cC7nV6zHIfyFP9EXjFa7YxKAVqPP4bxhhoLWW+z9JyCb6M/MREg59/RlmmXbmneIybB+YC/ay+yrffqEddDzwGvKxxDmzhc0tc80XVgblqFfgjwAAPubcGjAOl1wAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-undo {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4em8Dh0gAAAatJREFUWMPt1rFrFFEQBvDfGhACASshkL/ALpWVrSAKEQV5sIULWlgZNSgIFkGIVQ412gkBt1lYLERREFJqJRaW1oHAoZUQsDqwecWy7N3tbe6C4H2wxc682Zn3zTfvLXPM8b8j6RqYF+UCzsfnHBawGt3fMcAX7GEvS8NgKgXkRbmMxwg41TLsN0psZmnodyogL8pFPMIdLHUk7hA7eJKl4U/rAuKu3+HslFr/FZezNPSTFslX8QErDe4DvMVH/Iq9F7VwGpdwZUjsPtaSFjv/1vCBPjaxO0xcNbHejLpZrrlvJCMCT+JzA+2fcC1Lw+GE4l3CG1yIptfjCtiKoqtiJ0vD3aM0Py/K57iIMxgkQxat4EdN7e9xdRzlk+LEEPvDWvIDXJ928sYxjL36icWK+VaWhlezOIqbGFirJd/H7szugrwoX+D2BDEvszSsT5OBdfRaru/F9dPXQF6U27g/KnmWhgctxqyzBrZGMNGL/rHI0nDkKXiKexXTsywNGx0OnFbFNk3BRoWJXnw//j+ivCi32/S8CxPVNiWOAdUiJtXITIqYY45/Cn8B2D97FYW2H+IAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-wheel-pan {\n  background-image: url(\"data:image/png;base64,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\");\n}\n.bk-root .bk-tool-icon-wheel-zoom {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEgskILvMJQAABTtJREFUWMPdl1+MXVUVxn/fPvf2zrSFmUKnoBCUdjRoVaIxEpO2JhilMYBCtBQS2hejpg1Uo2NUrIFAoyGmtiE+GHwQGtvQJhqDmKYRBv+URFsFDNCSptH60DJTO3dKnX/33rM/H7rvsDu9M20fDMaVnGTvtb69z7fWXmvtc/TEzqd4OyXwNsv/FwFJQVI/sA14SZKRLOlPkr5TrVYXHz70quYkEEK4TtI2YAgYkrQthHDdhV5uuw+43/ZrwCbgRttgY/tjtrc0m83X3/f+D6ydnJhYcB4BSZcBA7aP2d4ELAGW2N5k+xgwkDB0IH19CGGH7R8B1aQeAf4KvAw0ku4K2zu7uru3ApdPEyiKohd4TNKjtjt5h6RHgccSNrddbvuHtm9Jqoak7xVF8WFgdavV+pSk5cCObNmXgK++85prCj3z28HKqZMnH7D9YAY4BvwujT8BvCuL1INX9vVt+dfwcCvNb7f9q2RuSfrGvWu/sL2Nf3LX7pzvj4ENSGBPVarVd4fRkZFltjdmoMGiKO4IIWwIIWwoiuIOYDDzeOPoyMiyFLkum7WJCMDztrcrTTrIRuAQZ6NcK1utL4dWq/VZoC8BhqvV6l1lWb4YYxyLMY6VZflitVq9CxhOmL60hhCKeYiV7WMKIXw9jT1HpXw3c+bOAKzOjJubzebJrKQCQLPZPClpc7bP6rMYKtjXth2OMf7tIkr11Wz8oQDc1Fb09vY+kQw1YAuwJY2nbUluAnCWpKkaFl6IQIzxivaR2SYA89sJVK/Xp2x32R6w/a30DNjuqtfrU0ArYecDCEqgLqm94T0dEm9mBG7PxkdDlkBnkhebgIezNQ8nHcCZPL9ijE1Jf/bZZoPtzbavmqNZLbf9tSxq+yoduuJ+SZ+zXSZyBXCqU+d8fvC5yRUrV+0G2j3g2hDCLyXd/+Su3QdnvP/zCuH72LWsgf2k0oHlH2c2odlkxcpVEdgr6aDtjyb8x20/J+mA7T9I6rL9SWA5dne2/GdXLl58qNJh398An85yTMA+4DOz8Dgu6Zu2dwJXJ91ltm8Gbp7Fgb+EEB4aHhpq5CEtACqVyr3AC0AlPS8k3TSmQ2YPhhBuS/1/LpmS9JTtNTHGfwBU2uUALARotVqniqJYH2Pck85pfavVaufAwnQvnHc0McaDKVptebN94QAnJB0EdtjekydyZXqjs/0ZgLIs/w6sy8bnYGYJ63pgERKC05JutT1kOwITwL9tvzlzUQUYB+Zjs2DBgu6xsbGJZHstByZbezregcBXeCsEz1bnzXt5anLyzLq71zDLxTRdVgemdx0fv2e2w5thO5DbiqL4oKT3ZKpnpyYnz+SY2ZpTAPZmJfdIrVZbNBNUq9UW2X4kU+2dcf53Aj1pj2PA7y/6m1DS00A9za9uNBq7iqJYBuoGdRdFsazRaOzKSqye1rTbaa/tlbYrqXQP2X4FIA9/J1l39xrC0v7+w5IeB8XkwS1lWe6TGJAYKMty31tfO4qSHl/a3384I3CDpI+kzC4lnRfrue6GytEjR8oQwlY73gC0L4qlth/q0M1/LYWtR48cKQF6enrC6dOnVwGLEpnxnp7en4+O1i/tszzGOCTpPmB7ahb57QUwBWyXdF+McWg6MScmuoA8OX8xOlpvXGz422XYTsB/SnpA0h7bX5R0WzI9HUL4qe2XbI+dk3xl+V7gxoztD5jRI+YK/zkEEokx2/uB/RdzIfUtueqVN04cXwF8G3iHY3z9Urw/j8ClyhsnjrcS2Vv/J/8NLxT+/zqBTkcxU/cfEkyEAu3kmjAAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-box-edit {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEg4QfHjM1QAAAGRJREFUWMNjXLhsJcNAAiaGAQYsDAwM/+lsJ+OgCwGsLqMB+D8o08CoA0YdMOqAUQewDFQdMBoFIyoN/B/U7YFRB7DQIc7xyo9GwbBMA4xDqhxgISH1klXbDYk0QOseEeOgDgEAIS0JQleje6IAAAAASUVORK5CYII=\");\n}\n.bk-root .bk-tool-icon-freehand-draw {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAADTElEQVRYCeWWTWwMYRjH/88721X1lZJIGxJxcEE4OOiBgzjXWh8TJKR76kWacOBGxdEJIdk4VChZI/phidRBHMRRIr7DSUiaSCRFRM3u88gz+o7Z6bBTdjmYZPf9eJ55fv/5zzvvDPC/H9QsA66Olo9Ga+/MdR+Ljm2/KQIULsz9FqItGdOfJKLhApLgVkiSCGODjWit7QpKWy+TNrFeXvzKVUT8NiTVaIgDcbiCFJ7GiT8WkARXAdYBK0Lbhi/CenArRNskuM7/tgNp4ArQ42dwjf3WY5gWTqC7O/NbNn2Xkfw/YwdSw/We14HP2IEZwX+y9cZ9SH0LmgFP7UCz4KkENBNeV0Cz4b8U8DfgKiDxMWwUXETqLvJpCQpXZfawbzS7t9v5pL19cHBwfja7YA0y/lyCM0+E5hv5+piZXwKYcF23as+37bTXsQVqgkL0p/34fHR7DcBtbetFsBmGDwMOJCggYG55yw7dMlk6DuC1Bdu2RsCU9TYWQq2IoGbsreZ5NzvEqfSBsIsIy8OTbcdgiRHeh4o8AFAEwDakbY2AaCCpH7V9aGhoUUUy3UyVbkPYFuYLDlUZH8XBpwxkK0Dbgxg5HcVi0ent7a0RULMIozaHBSMfF9b2SzdutFcFB2FkwMIJOG6qfteXOa1nHZ48tyefuwyfT9s6wtzZ3t7eZse2DR2I228TtHXzuWCx9g8MtK5cuHCZTH4tiHEOa4xFngvTyS8f35d6enomiCi4/foEXBkZaQuukChL4FYA2Whd7YcC4gEdW3CpdL3LtGAVCVYJywEyTpAuJKeMOKXZs/Bw947C50KhUFOG4cwz35cjWNBlHGeD53n3xsfHP/T19U1qciggar8Fa4I3PHobIotBWBtc2hSiChyZxVzM53Pv7FVH6Tp3uVy+g0r1ImD2GjIrQGYIxjnfuXTZGICS5k/bBwJoubwEFX4TLah9EXomJGMA3za+f9913Yl4TnzsDQ+vE6YTZOjHh4ngibstt1pzQwd04F0bPStEBpXqRoBeQ/AKghfBnOEKgS+Q7z91Xfdz/HGKg8Ox7z8iYD9z6wqTkZFgnvhMGP9VZ2or1XVkPM9z0mytSfVsHa1RLBZbLoyNzUnK+ydz3wC6I9x+lwbngwAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-poly-draw {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEjglo9eZgwAAAc5JREFUWMPt1zFrU1EUB/DfS4OmVTGDIChCP4BgnQXRxVHqIJUupp9AB8VBQcRBQUXIB9DWQoMRiXZzcnQSA34A7aAuHSJKkgo2LvfBrU3aJnlYkBy4vHcP557zP/9z3r33JdXa647N0kHSZd5Nn0rSxc8G3cXp85sMcnZZ8vge3osZ+l3vB8CWFA0iL14t79h210swAjACMAIwAjACkB90D/8/GchI9ve4nPwTBh5E9ws7OepzGWb9EddSn51Op9ZstadSg4VK1UKlKkmSDSMLALewiuNh/hVJq71Wxttmqz0dG88vPc+MgWP4grvYG3SLOBrZFFFrttqPe4HIDxh4GSei+98iSlusuYopXEAjBtEPA3tQwUpwluAbDm4TPJUz+BTW9l2Ce6G7L0X/Bw8D3T/7SKKIDzHg7QCcxjvcQAEtXAnrrg/RP0/DKPbqgcN4iVOR7gcO4dcQgRuoh7HSqwlP4n20m63jJu5n8MkWMYfP3UowhzdR8FU8w9iQwevBdyq3/27CMRzAE5yLuvsRLg+ZcR1nJ8YL81HWJUzGAPaFZwe/Q5MdyYDyNHgjzO90YyGHtVDncuiJchaHw8R4oREFV5qdiVmYLM3OgD9k5209/atmIAAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-point-draw {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEMEiERGWPELgAAA4RJREFUWMO1lr1uG1cQhb9ztdRSP7AF1QxgwKlcuZSqRC9gWUUUINWqTh5AnaFOnVPEteQmRuhCURqWsSqqc9IolREXdEvQBElxtdw7KURSFEVKu4w8wAKLxdw9Z+bMnRmZGXfZ29//II8th4WwGVNyIoQLYB5vxA9Caq04iUd9A+7ZlsNC2I7TdSd2hZXMJKlnTqp9jtl/GBaqoyQ0noFKpUIzBicYYc+DEFpxkglc4oVJa5gvDn8v1xV2irG3FM4NSVwjUKlUaMcpJhCGmSEJQ6QGD8M5WnHCd8+f3QCXpPLx8WNwv0j6Bm9FMK7FJ3WBE+R/2t7c/GBmFvSBrzRTCsyTDjXrxUgEMtpxynJYmJoBJ4VAybwVARgvL7Oik0okCodnKpVKX7P0leiVMb0VvbJT+upznK4vh0GIeQwwQStJkHQD3MwsCALTJRG7Qrdrj5m/djgYaIa0hlkRdJk26XEgC9txurccBtVW3IudBImmZuACUP+ZlIDBt9FKcubYNTcAH/X0RYM1E7utJPlqe+uZzPxUcEkiSS4sTT95n15Mud0xWC0o2PAWOCdK3KYZlFxfM+tHOcnMzNr1es18ug+cgsVjP4yBU/Ppfrter1m/+l0+zYygML1xRVHU7TSb1cSzBzoBzszsH+AMdJJ49jrNZjWKou6wBnwOzcyndBpNbuueURR1Dw8Pq35p9cc5p/Dy9Dypt7jXrtdGwQECS9NPhr6Gq6txUzNigE6zydLK6lTw12/KT4FGFEUfJX2YJNONq5tVs4ODA7sD/DnwJ/BoADZuE3tHFs12dna6d4C/BI6AlbyzI8ii2TTw12/KK33gb2cdXsNZoAntbZC2SeO4c9592k/5eNQbiwvFd1kJuFGwLJr1wSPg/SwpvyFBHufOeXcFeAlE97U/uCxOY+P3b+Bn4B3Q+L8EdJfD4a+/AbC4UBzPxiPg3wlHZquB28Cn2IuR9x3gr3uV4DbwfvSDOvi4uFA8BDZmIRHkjHpS9Ht9iRqd8+5G3g05mAGcQbsdiX5QJ428G7Kygo8XYdb1/K4NWVmjzkNge2sz84bs+ELmpDDLtqWsNZBXgvmw8CTtpWVMT7x5YWBjLARnwZfKQNYN2U2LPvrh+5nBt7c2M2/It9bArCTKR8eZN+SJ13AScPnoODeRdqNenH+wul5w2gUr2WUjMFAt8bZ/0axX/wNnv4H8vTFb1QAAAABJRU5ErkJggg==\");\n}\n.bk-root .bk-tool-icon-poly-edit {\n  background-image: url(\"data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gELFi46qJmxxAAABV9JREFUWMOdl19vFFUYxn9n9u9sCyylUIzWUoMQBAWCMdEEIt6xIRQSLIEKtvHe6AcA4yeQb7CAUNJy0daLeomJN8SEULAC2kBBapBKoLvbmdl/c14vdmY7u91tF95kknPOnHmf95znPc97Ro2OTeBbdjFDT3c32ZxVHUOE9kSMB0/m6ExuoJn1H+ur6Y+OTfD50SMN5168OgrAlyf7CfuD+z7+iDs3p8hkLUQ0iFQ/yFl5Nm/qonfHVva+s32Zw9GxCYILsZ08tpNfBhbs+1YN4OH9+7huGdECSBVfqUosbsllfmauBqiR+cCNwOr7AEo8pPHJnymXykhg5fUWjoQpl0vVvhZhbSzGoUOHqgBlt6B6uruj2Zy1E9jo0fhfeyL2x4Mnc8VErK0KUEOB64JSyptfG4RSytsJjUJVxw2lsFy3urL9nx1Qd25ObctkrVMi+jQivd7U2ZyV/3Hzpq7h3h1b/7p9Y0o8v8rwAbTWrGpSocN/FGDlbAI0Rl23PCBan0Ok158H9Ipwzi25A/Mzc9Gl/BYx/E4kYqC1NKRARNAaDCNUM27Z+Zr+ouXs0q4+LSLBHPYCFkTkC6uU39kwCdsS7WRKmaYUiAhdnZ3MPX2K4+QjQI+C94A93rMzm8ltMwyDeDzWjMZeEb2pYQDdW3vITU2jtUZ5QThOPgm8C7wP7J15OPsBsB3oWpGnVWisCeDS1VHj4vBI92+/3tgB7Ab2AruAXiDBK5oIOkhtkEYRNRuJhObrd8Dl9ewf4D5wG7hVLpen29vb5wzD+BrkbBMaL3d1dk5nsrnlFDTTFWAWmAZueWD3gCemGde2k2fw1Al1YXhEvjozoO49eczdqekrWmsc2zlrmvEKOGoW1GUjFLqSk2KpJrCLwyMCPAP+BO54QL8DM6YZX/ClsP9YnwKkXnIBP4jdIpJRpdJTCYdMwwi98KU0Hjc/dDILNyUcwTCWdOSMJ0TRmBktGRhLugu0xyLk7CIqVNm+0bGJptl1YXikD0grpY4Rjc4a8Fbgdab/6OGbAJeCUuyJnnHmZH9pbSyGuBXV8NUwlUpR1EWyixmSyTWEwqGlJ2Swbo2JXbAAfgDGgGQA9I1A9t1tlq0AxrXxn0ilUpw4fhQqYkH/sT41OTnJJwf2s6FjI5mshdYa7bqVR2uezr9MJmJt14FvGrh/O9D+e6UkM/xyCuCqEKCYnJyUTKFQrZDHjxzGshwWLQcRsOz8Hi85P23id0ug/XilAMLBmm4tPGdoaKjSH5+oAGrhwvBI9SjZTn4QSK9yenoD7dlrExPoJlXW8G8ytpNHxRKk02lGxsdRKFwXLNvx5yY94HQLGhGk4LFCYQSqaE0AwWM1eOoEbR0dKBSW7bC4mKuffxs4D/wCLKwQQPAUzIkslfp6cVomROWSolh0GjldAM4nzDi2k9/i5UAzC9aKfwNJ3zgJg9YEvN6+C7SHgKm69+sD7RfNnKTTaZRPQfAut4oFV//IS7gkcB34VlVo8kGzphlfB+DU+TfNGBpZtRastvrvARJmfMF28ge9sc2B9/PNnCilMIDwK6y8/ow/Ai4kvILTljAXvDvEvrqKSUs60KolzPjBxspavQD2tKqCAGF/Ba+xE/Wbilu54wZV8NEKF5fXzQHl/bh4hUsE0WAXSlDMYcQSrQXgCmsTseXHsJkNnjqBFGwKJaHsKlxtUHYVhbLCzr1kaOA4bcn1y1Swmb+iLpJKpVrfgdpfsiVVCYcgluwgnU7jEgJ4s5UkLFtWYyHyEg0/N1q1tmQH+YXnAMFr97Nmv3p+0QsHQRsF8qpBOE5+rb9Nkaj50tVQKjqh4OU3GNL/1/So3vuUgbAAAAAASUVORK5CYII=\");\n}\n");
export const bk_tool_icon_box_select = "bk-tool-icon-box-select";
export const bk_tool_icon_box_zoom = "bk-tool-icon-box-zoom";
export const bk_tool_icon_zoom_in = "bk-tool-icon-zoom-in";
export const bk_tool_icon_zoom_out = "bk-tool-icon-zoom-out";
export const bk_tool_icon_help = "bk-tool-icon-help";
export const bk_tool_icon_hover = "bk-tool-icon-hover";
export const bk_tool_icon_crosshair = "bk-tool-icon-crosshair";
export const bk_tool_icon_lasso_select = "bk-tool-icon-lasso-select";
export const bk_tool_icon_pan = "bk-tool-icon-pan";
export const bk_tool_icon_xpan = "bk-tool-icon-xpan";
export const bk_tool_icon_ypan = "bk-tool-icon-ypan";
export const bk_tool_icon_range = "bk-tool-icon-range";
export const bk_tool_icon_polygon_select = "bk-tool-icon-polygon-select";
export const bk_tool_icon_redo = "bk-tool-icon-redo";
export const bk_tool_icon_reset = "bk-tool-icon-reset";
export const bk_tool_icon_save = "bk-tool-icon-save";
export const bk_tool_icon_tap_select = "bk-tool-icon-tap-select";
export const bk_tool_icon_undo = "bk-tool-icon-undo";
export const bk_tool_icon_wheel_pan = "bk-tool-icon-wheel-pan";
export const bk_tool_icon_wheel_zoom = "bk-tool-icon-wheel-zoom";
export const bk_tool_icon_box_edit = "bk-tool-icon-box-edit";
export const bk_tool_icon_freehand_draw = "bk-tool-icon-freehand-draw";
export const bk_tool_icon_poly_draw = "bk-tool-icon-poly-draw";
export const bk_tool_icon_point_draw = "bk-tool-icon-point-draw";
export const bk_tool_icon_poly_edit = "bk-tool-icon-poly-edit";
//# sourceMappingURL=icons.js.map