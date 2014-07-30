file_exists() {
  if [[ -f "$1" ]]
  then
    return 0
    echo $1" exists!"
  else
    return 1
    echo $1" does not exist!"
  fi
}

directory_exists() {
  if [[ -d "$1" ]]
  then
    return 0
    echo $1" exists!"
  else
    return 1
    echo $1"  does not exist!"
  fi
}

load_config() {

  source $1

}
