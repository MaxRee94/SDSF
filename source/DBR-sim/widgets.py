"""All propeller widgets with a lower level of abstraction than the main window."""

from gui_utils import HContainer, VContainer, ListSettingWidget, KEYCODES
from PySide6 import QtCore, QtWidgets, QtGui


class SettingsContainer(QtWidgets.QWidget):
    """Contains publish settings widgets."""

    def __init__(self, model, parent=None):
        super(SettingsContainer, self).__init__(parent=parent)

        self.model = model
        self.setContentsMargins(4, 4, 4, 4)
        self.setPalette(QtGui.QPalette(QtGui.QColor(0, 0, 0)))

        self.layout = QtWidgets.QFormLayout(self)

        # Path widgets
        self.path_font = QtGui.QFont()
        self.path_widget = QtWidgets.QLabel("")
        self.path_widget.setFont(self.path_font)
        self.layout.addRow("Path:", self.path_widget)

        # Family widgets
        self.family_widget = FamilyContainer()
        self.layout.addRow("Family:", self.family_widget)

        # Source widgets
        source_container = QtWidgets.QWidget()
        source_layout = QtWidgets.QHBoxLayout(source_container)
        source_layout.setContentsMargins(0, 0, 0, 0)
        self.latest_button = QtWidgets.QPushButton("Latest")
        self.get_selected_button = QtWidgets.QPushButton("Selected")

        # TODO: Re-enable these buttons when their connections have been implemented
        self.latest_button.setEnabled(False)
        self.get_selected_button.setEnabled(False)

        self.source_line_edit = QtWidgets.QLineEdit()
        source_layout.addWidget(self.source_line_edit)
        source_layout.addWidget(self.latest_button)
        source_layout.addWidget(self.get_selected_button)
        self.layout.addRow("Source:", source_container)

        # Instance widgets
        instance_field = QtWidgets.QLineEdit()
        instance_field.setPlaceholderText("Enter instance name..")
        self.layout.addRow("Instance:", instance_field)

        # Instances widgets
        instances_container = InstancesContainer()
        self.layout.addRow("Instances:", instances_container)

        # Additional Families widgets
        self.families_listwidget = FamiliesContainer()
        self.layout.addRow("Additional\nfamilies:", self.families_listwidget)

        self.custom_instance_data_widget = QtWidgets.QTextEdit()
        self.custom_instance_data_widget.setPlaceholderText(
            "(Optional) Enter custom python dictionary..."
        )
        self.layout.addRow("Custom\nsettings:", self.custom_instance_data_widget)

        self.refresh()

    def refresh(self):
        print("refreshing..")


class SettingsWidget(QtWidgets.QWidget):
    """Displays and allows manipulation of publish settings."""

    def __init__(self, parent=None):
        super(SettingsWidget, self).__init__(parent=parent)

        layout = QtWidgets.QVBoxLayout(self)

        label = QtWidgets.QLabel("Publish settings")
        self.settings_container = SettingsContainer(None)

        layout.addWidget(label)
        layout.addWidget(self.settings_container)
        layout.addStretch()

    def set_file(self, file_id):
        self.data.set_file(file_id)

class FamilyContainer(QtWidgets.QComboBox):
    """Combobox widget to select a family from a list."""

    textChanged = QtCore.Signal()
    _PLACEHOLDER_TEXT = "Select.."

    def __init__(self, parent=None):
        super(FamilyContainer, self).__init__(parent=parent)

        self.currentTextChanged.connect(self._on_text_changed)

        self.populate()

    def populate(self):
        """Populate the list of family options."""
        pass

    def text(self):
        """Return the text associated with the currently selected item."""
        return self.currentText()

    def setText(self, text):
        self.setCurrentText(text)

    def _on_text_changed(self):
        self.textChanged.emit()

    def setPlaceholderText(self, _):
        self.setCurrentIndex(self.count() - 1)


class FamiliesContainer(ListSettingWidget):
    """Widget for selecting multiple additional families."""

    _NEW_ITEM_FIELD = FamilyContainer
    _PLACEHOLDER_TEXT = "Select.."

    def __init__(self, parent=None):
        super(FamiliesContainer, self).__init__(parent=parent)

        self._IS_ITEM_NAME_VALID = self.is_item_name_valid

    def is_item_name_valid(self, name):
        for i in range(self.listwidget.count()):
            item = self.listwidget.item(i)
            item_txt = item.text()
            if item_txt == name:
                return False
        else:
            if name == self._PLACEHOLDER_TEXT:
                return False

        return True

    def refresh(self):
        super(FamiliesContainer, self).refresh()

        if self.new_item_field.text() == self._PLACEHOLDER_TEXT:
            self.add_btn.setEnabled(False)
        else:
            self.add_btn.setEnabled(True)

    def get_items_texts(self):
        items = []
        for i in range(self.listwidget.count()):
            item = self.listwidget.item(i)
            item_txt = item.text()
            items.append(item_txt)

        return items


class InstancesContainer(ListSettingWidget):
    """Widget for selecting multiple instances."""

    _PLACEHOLDER_TEXT = "Enter instance name.."


