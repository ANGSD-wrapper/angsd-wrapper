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

file_not_empty() {
  if [[ -s "$1"]]
  then
    return 0
    echo $1" not empty!"
  else
    return 1
    echo $1" empty!"
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

variable_exists() {
  if [[ -z "$1"]]
  then
    return 0
    echo $1" exists!"
  else
    return 1
    echo $1" does not exist!"
  fi
}

load_config() {

  source $1

}
