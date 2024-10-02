storm_path=../storm/build/bin/storm-conv
model_path=../mdp-models/action-models-prism
jani_path=../mdp-models/action-models-jani
for m in $model_path/pnueli-zuck*.prism;
do
	$storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/pnueli-zuck.props --globalvars
done

for m in $model_path/consensus.*.prism;
do
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/consensus.props --globalvars
done

for m in $model_path/csma.*.prism;
do
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/csma.props --globalvars
done

for m in $model_path/firewire_dl.prism;
do
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/firewire_dl.props --globalvars
done


for m in $model_path/mer.prism;
do
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/mer.props --globalvars
done

for m in $model_path/pacman.prism;
do
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/pacman.props --globalvars
done

for m in $model_path/philosophers.*.prism;
do
  props=$model_path/$(basename $m).props
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop "${props/.prism.props/.props}" --globalvars
done

for m in $model_path/zeroconf_dl.prism;
do
  $storm_path --prism $m --tojani $jani_path/$(basename $m).jani --prop $model_path/zeroconf_dl.props --globalvars
done

for file in $jani_path/*.prism.jani
do
	mv "$file" "${file%.prism.jani}.jani"
done

