check_file() {
  if [ -f "$1" ]
  then
    return 0
  else
    return 1
  fi
}

check_directory()
  if [ -d "$1" ]
  then
    return 0
  else
    return 1
  fi
}